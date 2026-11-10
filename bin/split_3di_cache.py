#!/usr/bin/env python3
"""
Partition a protein FASTA by whether a previously-computed ProstT5 3Di
prediction already exists in the persistent cache, so ProstT5 (a multi-hour
GPU step) is never re-run for a protein it's already predicted -- on any
species, any project, since TransDecoder headers (gene_id~~transcript_id,
with ORF suffix) are globally unique.

Cache layout: <cache_dir>/<accession>.3di, one file per protein, holding
just the predicted 3Di-alphabet sequence as plain text with no FASTA
header -- the same layout `foldseek unpackdb --unpack-name-mode 1
--unpack-suffix .3di` produces, where <accession> is exactly the FASTA
header up to the first whitespace (foldseek's own key convention, matching
what it records in a DB's .lookup file).

Usage:
    split_3di_cache.py <pep.fasta> <cache_dir> \\
        <cached_aa.fasta> <cached_3di.fasta> <needs_prostt5.fasta>

Writes cached_aa.fasta and cached_3di.fasta with matching entries in the
same order (so a plain `createdb` on each -- no --prostt5-model -- produces
two DBs whose internal ids line up 1:1, letting the 3di one be relabelled
as a "_ss" object without ever invoking ProstT5), and needs_prostt5.fasta
with everything not yet cached.
"""

import sys
import os


def read_fasta_records(path):
    """Yield (header, [sequence lines]) for each record, in order."""
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
    if len(sys.argv) != 6:
        sys.exit(
            f"Usage: {sys.argv[0]} <pep.fasta> <cache_dir> "
            "<cached_aa.fasta> <cached_3di.fasta> <needs_prostt5.fasta>"
        )

    pep_fasta, cache_dir, cached_aa_path, cached_3di_path, needs_path = sys.argv[1:6]

    n_cached = 0
    n_needs = 0
    with open(cached_aa_path, "w") as cached_aa, \
         open(cached_3di_path, "w") as cached_3di, \
         open(needs_path, "w") as needs:
        for header, seq_lines in read_fasta_records(pep_fasta):
            accession = header[1:].split()[0]
            cache_path = os.path.join(cache_dir, f"{accession}.3di")
            if os.path.isfile(cache_path):
                with open(cache_path) as fh:
                    threedi_seq = fh.read().strip()
                cached_aa.write(header + "\n")
                cached_aa.write("\n".join(seq_lines) + "\n")
                cached_3di.write(header + "\n")
                cached_3di.write(threedi_seq + "\n")
                n_cached += 1
            else:
                needs.write(header + "\n")
                needs.write("\n".join(seq_lines) + "\n")
                n_needs += 1

    sys.stderr.write(
        f"split_3di_cache: {n_cached} proteins reused from cache, "
        f"{n_needs} need fresh ProstT5 prediction\n"
    )


if __name__ == "__main__":
    main()
