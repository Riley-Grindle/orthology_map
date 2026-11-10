#!/usr/bin/env python3
"""
Expand i-ADHoRe's multiplicons.txt (one row per syntenic *block*) into a
gene-level table: for every (non-redundant) block's segment, every gene
physically between the segment's first and last gene (inclusive, in the
real per-chromosome gene-order file we built and fed to i-ADHoRe) gets a
row, tagged with that block.

Unlike i-ADHoRe's own list_elements.txt, this includes genes that AREN'T
part of the matched synteny signal -- "gap" genes sitting between two
actual chain members with no homolog of their own in the other genome.
Each row's `contributes` column says which is which: 1 if the gene is one
of the actual homology/anchor-chain members i-ADHoRe reported in
list_elements.txt for this segment, 0 if it's just physically in between.

i-ADHoRe's own `position` numbers in list_elements.txt don't correspond to
simple file/line order (confirmed empirically against the real .lst files
-- the gaps between consecutive reported positions don't match the gaps
between the corresponding genes' real line numbers, likely due to internal
tandem-duplicate handling that isn't exposed in any output file), so they
can't be used to reconstruct "everything in between". Instead, `position`
here is this script's OWN 0-indexed offset from the segment's first gene,
based on the real gene order in the .lst files -- the same files used to
build i-ADHoRe's input, so this order is unambiguous.

Each row's genome/list are resolved via segments.txt (uniform for nested/
child blocks too). A gene can legitimately appear in more than one block's
output rows if blocks overlap -- that's not deduplicated here, since a
gene genuinely can belong to multiple syntenic contexts.

parent_block_id (empty for a top-level block) and depth (0 for a top-level
block, 1 for its direct children, 2 for their children, etc., walking
multiplicons.txt's `parent` pointers) are also carried per row. This is
NOT the same as multiplicons.txt's own `level` column -- `level` is
i-ADHoRe's internal N-way multiplicon degree (2 for an ordinary pairwise
block, 3+ for one involving an extra paralogous relationship), not a
nesting counter, and doesn't start at 0.

Usage:
    multiplicon_to_genes.py <i-adhore output dir> <output.tsv> \\
        <query_genome_id> <query_genelists_dir> \\
        <ref_genome_id> <ref_genelists_dir>

Expects <output dir> to contain multiplicons.txt, segments.txt, and
list_elements.txt (i-ADHoRe's standard output set). <query_genelists_dir>/
<ref_genelists_dir> are directories of per-chromosome .lst files (one gene
per line, in genomic order, id+strand with no separator), keyed by
filename <list>.lst -- the same files i-ADHoRe itself was given.
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


def load_gene_order(lst_path):
    """Return (ordered [(gene_id, strand), ...], {gene_id: index}) for one .lst file."""
    genes = []
    with open(lst_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            # entries are "<gene_id~~transcript_id><+|-|.>" with no separator
            genes.append((line[:-1], line[-1]))
    index = {gene_id: i for i, (gene_id, _strand) in enumerate(genes)}
    return genes, index


def main():
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <i-adhore output dir> <output.tsv> "
            "<query_genome_id> <query_genelists_dir> <ref_genome_id> <ref_genelists_dir>"
        )

    out_dir, out_path, query_id, query_dir, ref_id, ref_dir = sys.argv[1:7]
    multiplicons_path = os.path.join(out_dir, "multiplicons.txt")
    segments_path = os.path.join(out_dir, "segments.txt")
    list_elements_path = os.path.join(out_dir, "list_elements.txt")

    genelists_dir = {query_id: query_dir, ref_id: ref_dir}

    # multiplicons.txt columns: id genome_x list_x parent genome_y list_y
    # level number_of_anchorpoints profile_length begin_x end_x begin_y
    # end_y is_redundant
    valid_blocks = set()
    parent_of = {}
    for fields in read_tsv(multiplicons_path):
        block_id, parent_id, is_redundant = fields[0], fields[3], fields[13]
        # Keep every row's parent pointer, even a redundant one, so depth
        # can still be walked correctly if a block's ancestor happens to be
        # one that's excluded from the output as redundant.
        parent_of[block_id] = parent_id
        if is_redundant != "1":
            valid_blocks.add(block_id)

    depth_cache = {}

    def depth_of(block_id):
        """0 for a top-level block, 1 for its direct children, etc."""
        if block_id in depth_cache:
            return depth_cache[block_id]
        depth, current, seen = 0, block_id, {block_id}
        while parent_of.get(current):
            current = parent_of[current]
            if current in seen:
                break  # guard against malformed/cyclic parent data
            seen.add(current)
            depth += 1
        depth_cache[block_id] = depth
        return depth

    # segments.txt columns: id multiplicon genome list first last order
    # Map each segment id -> (block id, genome, list, first gene, last gene),
    # restricted to blocks kept above.
    segment_info = {}
    for fields in read_tsv(segments_path):
        segment_id, block_id, genome, chrom_list, first_gene, last_gene = fields[0:6]
        if block_id in valid_blocks:
            segment_info[segment_id] = (block_id, genome, chrom_list, first_gene, last_gene)

    # list_elements.txt columns: id segment gene position orientation
    # Per segment, which genes i-ADHoRe reported as actual chain members
    # (anything else between first and last is a "gap" gene, contributes=0).
    contributing = {}
    for fields in read_tsv(list_elements_path):
        segment_id, gene = fields[1], fields[2]
        contributing.setdefault(segment_id, set()).add(gene)

    gene_order_cache = {}

    def gene_order_for(genome, chrom_list):
        key = (genome, chrom_list)
        if key not in gene_order_cache:
            lst_path = os.path.join(genelists_dir[genome], f"{chrom_list}.lst")
            gene_order_cache[key] = load_gene_order(lst_path)
        return gene_order_cache[key]

    with open(out_path, "w") as out:
        out.write(
            "block_id\tparent_block_id\tdepth\tgenome\tlist\tgene\t"
            "position\torientation\tcontributes\n"
        )
        for segment_id, (block_id, genome, chrom_list, first_gene, last_gene) in segment_info.items():
            genes, index = gene_order_for(genome, chrom_list)
            if first_gene not in index or last_gene not in index:
                continue
            start, end = index[first_gene], index[last_gene]
            if start > end:
                start, end = end, start
            members = contributing.get(segment_id, set())
            parent_id = parent_of.get(block_id, "")
            depth = depth_of(block_id)
            for offset, i in enumerate(range(start, end + 1)):
                gene_id, strand = genes[i]
                contributes = 1 if gene_id in members else 0
                out.write(
                    f"{block_id}\t{parent_id}\t{depth}\t{genome}\t{chrom_list}\t"
                    f"{gene_id}\t{offset}\t{strand}\t{contributes}\n"
                )


if __name__ == "__main__":
    main()
