#!/usr/bin/env python3
"""
Build a minimal foldseek/MMseqs2-format sequence DB directly from a FASTA
file, without going through foldseek's own `createdb` -- which only
accepts plain sequence/text input via --prostt5-model (triggering a real
ProstT5 prediction) and otherwise always takes the structure-parsing path,
even on a FASTA with no PDB/mmCIF content at all (confirmed empirically:
it exits non-zero with "No structures found in given input" and writes
empty output files).

This is needed to reuse a previously-cached ProstT5 prediction (or any
other pre-existing text, treated as opaque sequence content) as a foldseek
DB object -- e.g. to relabel as a "_ss" (3Di) object -- without invoking
ProstT5 again. The on-disk format doesn't care what alphabet the text is
in, so this works for either amino-acid or 3Di-alphabet FASTA input.

Format (reverse-engineered from real foldseek output, byte-for-byte):
  <out_prefix>          -- data file: for each entry in input order,
                           "<sequence text>\n\0", concatenated with no
                           further separator
  <out_prefix>.index    -- "<0-indexed key>\t<byte offset>\t<entry length,
                           including the trailing \n\0>\n" per entry
  <out_prefix>.dbtype   -- 4 bytes: 0x00 0x00 0x00 0x00 (sequence-type,
                           the same code foldseek uses for both AA and 3Di
                           sequence objects)
  <out_prefix>.lookup   -- "<key>\t<accession>\t0\n" per entry, where
                           accession is the FASTA header up to the first
                           whitespace -- foldseek's own key convention
  <out_prefix>.source   -- "0\t<fasta path>\n" (single line; every entry
                           comes from the one input file, source index 0)
  <out_prefix>_h(.index/.dbtype) -- the header-DB companion object, same
                           format as the main data file but holding each
                           entry's full header text (no leading ">")
                           instead of its sequence -- always produced
                           alongside a sequence DB, the same as createdb's
                           own output.

Usage:
    fasta_to_seqdb.py <input.fasta> <out_prefix>
"""

import sys


def read_fasta_records(path):
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            elif header is not None:
                seq_lines.append(line)
    if header is not None:
        yield header, seq_lines


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.fasta> <out_prefix>")

    fasta_path, out_prefix = sys.argv[1:3]

    with open(out_prefix, "wb") as data_fh, \
         open(f"{out_prefix}.index", "w") as index_fh, \
         open(f"{out_prefix}.lookup", "w") as lookup_fh, \
         open(f"{out_prefix}_h", "wb") as h_data_fh, \
         open(f"{out_prefix}_h.index", "w") as h_index_fh:
        offset = 0
        h_offset = 0
        key = 0
        for header, seq_lines in read_fasta_records(fasta_path):
            header_text = header[1:]
            accession = header_text.split()[0]
            seq = "".join(seq_lines)

            entry = (seq + "\n").encode() + b"\0"
            data_fh.write(entry)
            index_fh.write(f"{key}\t{offset}\t{len(entry)}\n")
            lookup_fh.write(f"{key}\t{accession}\t0\n")
            offset += len(entry)

            h_entry = (header_text + "\n").encode() + b"\0"
            h_data_fh.write(h_entry)
            h_index_fh.write(f"{key}\t{h_offset}\t{len(h_entry)}\n")
            h_offset += len(h_entry)

            key += 1

    for suffix in ("", "_h"):
        with open(f"{out_prefix}{suffix}.dbtype", "wb") as fh:
            fh.write(b"\x00\x00\x00\x00")

    with open(f"{out_prefix}.source", "w") as fh:
        fh.write(f"0\t{fasta_path}\n")


if __name__ == "__main__":
    main()
