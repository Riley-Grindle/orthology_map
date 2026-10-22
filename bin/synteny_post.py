#!/usr/bin/env python3
"""
Convert i-ADHoRe's synteny_std.tsv (QUERY_ID, REF_ID, SCORE) into the
<QUERY, MATCH, GENE, SPECIES> format tbl_2_json.py requires for every other
1v1/1vM evidence source. Mirrors the Ensembl-id lookup that tree_post.R's
pathway uses (ensembl_id_2_gene_symbl.py), reusing the same supplementary
data (prefix_2_file.json / prefix_2_species.json baked into the post_proc
image at /sup_data/, and the per-species gene-symbol JSON files under
--ensembl_data) so REF_ID gets resolved to a real gene symbol + species name.

REF_ID is "gene_id~~transcript_id" (TransDecoder's --gene_trans_map
convention); only the gene_id half is looked up.

Usage:
    synteny_post.py <synteny_std.tsv> <ensembl_data_dir> <prefix_2_file.json> \
                     <prefix_2_species.json> <output.tsv>
"""

import sys
import json
import os


def get_prefix(ensembl_id):
    """Same convention as ensembl_id_2_gene_symbl.py's get_prefix()."""
    if ensembl_id.startswith("MGP_"):
        end = ensembl_id.index("_", 4)
        return ensembl_id[: end + 1]

    alpha_run = 0
    while alpha_run < len(ensembl_id) and ensembl_id[alpha_run].isalpha():
        alpha_run += 1

    if alpha_run < 2:
        return None

    return ensembl_id[: alpha_run - 1]


def main():
    synteny_tsv, ensembl_data, prefix_2_file, prefix_2_species, out_path = sys.argv[1:6]

    with open(prefix_2_file) as fh:
        prefix_dict = json.load(fh)
    with open(prefix_2_species) as fh:
        species_dict = json.load(fh)

    gene_maps = {}

    def lookup(gene_id):
        prefix = get_prefix(gene_id)
        if prefix is None or prefix not in prefix_dict:
            return "NA", "NA"

        if prefix not in gene_maps:
            data_file = os.path.join(ensembl_data, prefix_dict[prefix])
            with open(data_file) as fh:
                gene_maps[prefix] = json.load(fh)

        symbol = gene_maps[prefix].get(gene_id, "NA")
        species = species_dict.get(prefix, "NA").strip()
        return symbol, species

    with open(synteny_tsv) as fh, open(out_path, "w") as out:
        out.write("QUERY\tMATCH\tGENE\tSPECIES\n")
        next(fh, None)  # header
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            query_id, ref_id = fields[0], fields[1]
            ref_gene = ref_id.split("~~")[0]
            symbol, species = lookup(ref_gene)
            out.write(f"{query_id}\t{ref_id}\t{symbol}\t{species}\n")


if __name__ == "__main__":
    main()
