#!/usr/bin/env python3
"""
Convert a tool's standardized <QUERY, MATCH, GENE, SPECIES> TSV into a
{query_gene: {gene_symbol: [species, ...]}} JSON.

Continuation rows: when a query gene has multiple matches, the post-
processing R scripts (ortho_f_post.R / ortho_l_post.R / tree_post.R /
eggnog_post.R) only write the QUERY id on the first row of that gene's
block; later rows leave QUERY blank to mean "same gene as the row above".
This is tracked via `key` across rows instead of re-deriving it each time.
"""

import sys
import json


def parse_hits(tbl_lst):
    genes = {}
    key = None

    for lineno, line in enumerate(tbl_lst[1:], start=2):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            sys.exit(
                f"ERROR line {lineno}: expected >=4 tab-separated fields, "
                f"got {len(fields)}: {line!r}"
            )

        query, _match, symbol, species = fields[0], fields[1], fields[2], fields[3]

        if "~" in query:
            key = query.split("~")[0]
            genes.setdefault(key, [])
        elif key is None:
            sys.exit(
                f"ERROR line {lineno}: QUERY field has no gene delimiter "
                "('~') and no earlier row established a gene key -- the "
                f"first data row of a std.tsv must carry a real query id: {line!r}"
            )

        if symbol not in ("NA", "NO GENE SYMBOL"):
            genes[key].append((symbol.upper(), species.strip()))

    return genes


def join_hits(match_list):
    tmp_dict = {}
    for symbol, species in match_list:
        tmp_dict.setdefault(symbol, set()).add(species)
    return {symbol: list(species_set) for symbol, species_set in tmp_dict.items()}


def main():
    with open(sys.argv[1]) as fh:
        file_lst = fh.readlines()

    genes = parse_hits(file_lst)
    for gene, matches in genes.items():
        genes[gene] = join_hits(matches)

    if len(sys.argv) > 3:
        file_name = f"{sys.argv[2]}.{sys.argv[3]}.query_2_matches.json"
        genes["tool"] = sys.argv[3]
    else:
        file_name = f"{sys.argv[2]}.query_2_matches.json"
        genes["tool"] = sys.argv[2]

    with open(file_name, "w") as output:
        json.dump(genes, output, indent=4)


if __name__ == "__main__":
    main()
