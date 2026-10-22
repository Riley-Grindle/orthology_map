#!/usr/bin/env python3
"""
Partition a protein FASTA into two sets based on whether the user already
supplied a PDB file for each protein, for structural-comparison "alphafold"
mode: proteins with a pre-folded structure get symlinked into a directory
for foldseek to ingest as real structures; everything else gets written to
a fallback FASTA for foldseek to run ProstT5 3Di-prediction on instead.

FASTA headers follow this pipeline's TransDecoder convention:
    >gene_id~~transcript_id.pN  (optional ORF suffix after the transcript id)
The transcript id (without the ORF suffix) is what's expected to name each
supplied PDB file: <transcript_id>.pdb

Usage:
    split_pdb_coverage.py <pep.fasta> <pdb_dir> <has_pdb_outdir> <missing.fasta>
"""

import os
import re
import sys


def transcript_id_from_header(header):
    """Extract the transcript id from a TransDecoder-style FASTA header."""
    name = header[1:].split()[0]
    if "~~" in name:
        name = name.split("~~", 1)[1]
    return re.sub(r"\.p\d+$", "", name)


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
    if len(sys.argv) != 5:
        sys.exit(f"Usage: {sys.argv[0]} <pep.fasta> <pdb_dir> <has_pdb_outdir> <missing.fasta>")

    pep_fasta, pdb_dir, has_pdb_outdir, missing_fasta = sys.argv[1:5]
    os.makedirs(has_pdb_outdir, exist_ok=True)

    n_has_pdb = 0
    n_missing = 0
    with open(missing_fasta, "w") as missing_fh:
        for header, seq_lines in read_fasta_records(pep_fasta):
            transcript_id = transcript_id_from_header(header)
            pdb_path = os.path.join(pdb_dir, f"{transcript_id}.pdb")
            if os.path.isfile(pdb_path):
                link_path = os.path.join(has_pdb_outdir, f"{transcript_id}.pdb")
                if not os.path.exists(link_path):
                    os.symlink(os.path.abspath(pdb_path), link_path)
                n_has_pdb += 1
            else:
                missing_fh.write(header + "\n")
                missing_fh.write("\n".join(seq_lines) + "\n")
                n_missing += 1

    sys.stderr.write(
        f"split_pdb_coverage: {n_has_pdb} proteins have a supplied PDB, "
        f"{n_missing} fall back to ProstT5\n"
    )


if __name__ == "__main__":
    main()
