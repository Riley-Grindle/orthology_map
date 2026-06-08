#!/usr/bin/env python3
"""
Build i-ADHoRe gene-list files from a GTF.

One .lst file is written per chromosome/scaffold.  Each file lists genes in
genomic order with their strand appended directly to the ID (no space):

    TRINITY_DN1_c0_g1_i1+
    ENST00000380152.7-

IDs are transcript_id values so they match the headers in TransDecoder's
final.transdecoder.pep (which also uses transcript_id after stripping .pN).
Genes are de-duplicated and sorted by start coordinate within each scaffold.

Usage:
    gtf_to_genelists.py <input.gtf> <output_dir>
"""

import sys
import os
import re
from collections import defaultdict


def parse_attrs(attr_str):
    attrs = {}
    for m in re.finditer(r'(\w+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def gtf_to_genelists(gtf_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # chrom -> {tx_id: (start, strand)}
    chroms = defaultdict(dict)

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue

            feature = fields[2]
            chrom   = fields[0]
            start   = int(fields[3])
            strand  = fields[6]
            attrs   = parse_attrs(fields[8])

            if feature in ('transcript', 'gene'):
                tx_id = attrs.get('transcript_id') or attrs.get('gene_id', '')
            elif feature == 'exon':
                tx_id = attrs.get('transcript_id', '')
            else:
                continue

            if tx_id and tx_id not in chroms[chrom]:
                chroms[chrom][tx_id] = (start, strand)

    if not chroms:
        sys.exit(f"ERROR: no transcript/exon features found in {gtf_path}")

    for chrom, genes in sorted(chroms.items()):
        sorted_genes = sorted(genes.items(), key=lambda kv: kv[1][0])
        lst_path = os.path.join(out_dir, f"{chrom}.lst")
        with open(lst_path, 'w') as fh:
            for tx_id, (_, strand) in sorted_genes:
                fh.write(f"{tx_id}{strand}\n")

    n = sum(len(v) for v in chroms.values())
    print(f"Wrote {n} genes across {len(chroms)} scaffolds to {out_dir}/",
          file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.gtf> <output_dir>")
    gtf_to_genelists(sys.argv[1], sys.argv[2])
