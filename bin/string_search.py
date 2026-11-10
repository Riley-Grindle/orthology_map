#!/usr/bin/env python3

# Python3 code to implement iterative Binary
# Search.

import sys
import json


# Driver Code
if __name__ == '__main__':
    name = sys.argv[1]
    file = open(name, "r")
    dict_ids = json.load(file)


    target_name = sys.argv[2]
    file_target = open(target_name, "r")
    target_list = file_target.readlines()


    # Dictionary Search
    #
    # NOTE: dict_ids[key] used to be a bare lookup with no fallback, so any
    # STRING-DB id not covered by the local string_db.json subset (e.g. a
    # less common organism's protein) crashed the whole run with a
    # KeyError instead of just leaving that one row unresolved. Match the
    # "NA" fallback convention used elsewhere in this pipeline (see
    # ensembl_id_2_gene_symbl.py) so one missing id doesn't take down the
    # entire run.
    gene_symbols = list()
    species_list = list()
    for i in range(len(target_list)):

        if i > 0:
            key = target_list[i].strip()
            try:
                result = dict_ids[key].strip()
            except KeyError:
                result = "NA"
            gene_symbols.append(result)
            tmp_lst = [*target_list[i]]
            try:
                index = tmp_lst.index(".")
                species_list.append(target_list[i][:index])
            except ValueError:
                species_list.append("NA")

    # NOTE: this used to build long_string/spec_long_string via repeated
    # string concatenation (`x = x + ...`) in a loop -- the same O(n^2)
    # pattern fixed in ensembl_id_2_gene_symbl.py. Appending to a list and
    # joining once at the end is O(n).
    gene_lines = ["GENE"]
    species_lines = ["SPECIES"]
    for i in range(len(gene_symbols)):
        gene_lines.append(gene_symbols[i])
        species_lines.append(species_list[i])

    with open("gene_symbols.txt", "w") as gene_file:
        gene_file.write("\n".join(gene_lines) + "\n")

    with open("taxa_ids.txt", "w") as species_file:
        species_file.write("\n".join(species_lines) + "\n")
