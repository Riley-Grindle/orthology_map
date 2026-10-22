#!/usr/bin/env python3
"""
Expand i-ADHoRe's multiplicons.txt (one row per syntenic *block*) into a
gene-level table: one row per gene that falls inside any (non-redundant)
block's segment, tagged with which block it belongs to.

multiplicons.txt only gives each block's boundary as a gene-order range
(begin_x/end_x, begin_y/end_y) per genome, and for nested/child blocks
(those with a non-empty `parent`, from i-ADHoRe's generic N-way profiling --
even a 2-genome comparison can produce these as paralogous regions that
share synteny with the same reference locus) one side's genome/list isn't
given directly in that row at all. segments.txt resolves this uniformly: it
has a row per block per genome, with genome/list and the actual first/last
gene IDs bounding that segment, regardless of nesting. list_elements.txt
then gives every individual gene within a segment (not just the two
boundary genes), which is what makes this expansion "all genes in the
block" rather than just its two endpoints.

So the join is: multiplicons.txt (which block ids to keep) -> segments.txt
(block id -> its segment ids, with genome/list) -> list_elements.txt
(segment id -> every gene in it). A gene can legitimately appear in more
than one block's output rows if blocks overlap -- that's not deduplicated
here, since a gene genuinely can belong to multiple syntenic contexts.

Usage:
    multiplicon_to_genes.py <i-adhore output dir> <output.tsv>

Expects <output dir> to contain multiplicons.txt, segments.txt, and
list_elements.txt (i-ADHoRe's standard output set).
"""

import sys
import os


def read_tsv(path):
    """Yield each data row (header skipped) as a list of fields."""
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            line = line.rstrip("\n")
            if line:
                yield line.split("\t")


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <i-adhore output dir> <output.tsv>")

    out_dir, out_path = sys.argv[1], sys.argv[2]
    multiplicons_path = os.path.join(out_dir, "multiplicons.txt")
    segments_path = os.path.join(out_dir, "segments.txt")
    list_elements_path = os.path.join(out_dir, "list_elements.txt")

    # multiplicons.txt columns: id genome_x list_x parent genome_y list_y
    # level number_of_anchorpoints profile_length begin_x end_x begin_y
    # end_y is_redundant
    valid_blocks = set()
    for fields in read_tsv(multiplicons_path):
        block_id, is_redundant = fields[0], fields[13]
        if is_redundant != "1":
            valid_blocks.add(block_id)

    # segments.txt columns: id multiplicon genome list first last order
    # Map each segment id -> (block id, genome, list), restricted to blocks
    # kept above.
    segment_info = {}
    for fields in read_tsv(segments_path):
        segment_id, block_id, genome, chrom_list = fields[0], fields[1], fields[2], fields[3]
        if block_id in valid_blocks:
            segment_info[segment_id] = (block_id, genome, chrom_list)

    # list_elements.txt columns: id segment gene position orientation
    with open(list_elements_path) as fh, open(out_path, "w") as out:
        out.write("block_id\tgenome\tlist\tgene\tposition\torientation\n")
        next(fh, None)  # header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            _id, segment_id, gene, position, orientation = line.split("\t")
            info = segment_info.get(segment_id)
            if info is None:
                continue
            block_id, genome, chrom_list = info
            out.write(f"{block_id}\t{genome}\t{chrom_list}\t{gene}\t{position}\t{orientation}\n")


if __name__ == "__main__":
    main()
