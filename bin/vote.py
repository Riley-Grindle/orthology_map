#!/usr/bin/env python3
"""
Union every tool's <gene>.query_2_matches.json into voted_orthologs.json,
tagging each candidate match with the (species, tool) pairs that support it.

This performs no filtering -- every match any tool proposed is kept. The
inclusion/confidence decision (best_match, confidence tier) is made
afterwards by rank_matches.py, which reads the fully unioned evidence this
script produces.
"""

import sys
import json


def extract_gene_id(query_list):
    gene_dict = {}
    for line in query_list:
        gene_dict[line[1:line.index("~")]] = {}
    return gene_dict


def combine(agg, comp):
    for key, value in comp.items():
        if key == "tool":
            continue
        if key not in agg:
            agg[key] = {}
        for match, tuples in value.items():
            if match in agg[key]:
                agg[key][match] = agg[key][match] + tuples
            else:
                agg[key][match] = tuples
    return agg


def tuplify(comp, tool):
    for key, value in comp.items():
        if key == "tool":
            continue
        for gene, info in value.items():
            comp[key][gene] = [tuple([taxa, tool]) for taxa in info]
    return comp


def main():
    with open(sys.argv[1], "r") as query_file:
        query_list = query_file.readlines()

    q_genes = extract_gene_id(query_list)

    with open(sys.argv[3], "r") as comparing_file:
        comparison = json.load(comparing_file)

    with open(sys.argv[2], "r") as aggregate_file:
        try:
            aggregate = json.load(aggregate_file)
        except json.decoder.JSONDecodeError:
            aggregate = q_genes

    formatted_comp = tuplify(comparison, comparison["tool"])
    new_aggregate = combine(aggregate, formatted_comp)

    with open(sys.argv[2], "w") as outfile:
        json.dump(new_aggregate, outfile, indent=4)


if __name__ == "__main__":
    main()
