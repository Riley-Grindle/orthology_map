#!/usr/bin/env python3
"""
Replace tree_std.csv's QUERY ids with the true original ids, via an exact
lookup built by build_treegrafter_id_map.py.

tree_post.R recovers TreeGrafter's sanitized query id and tries to restore
the original "GENE~~transcript.pN" delimiters with gsub("_", "|", ...) --
converting every underscore to "|" -- on the assumption that every
underscore in the sanitized id came from a "~~" or "." that TreeGrafter
itself sanitized away. That assumption breaks whenever the real gene or
transcript id contains a natural underscore, silently corrupting the QUERY
id and breaking the join with the other tools' evidence for that gene. This
script undoes that specific substitution (turning "|" back into "_") to
recover TreeGrafter's sanitized id exactly, then looks it up directly
instead of guessing where the real delimiters were.

Continuation rows (QUERY blank -- tree_post.R leaves it blank/whitespace to
mean "same query gene as the row above") are passed through unchanged.
"""

import csv
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <tree_std.csv> <id_map.tsv> <output.csv>")

    csv_path, map_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    id_map = {}
    with open(map_path) as fh:
        for line in fh:
            sanitized, original = line.rstrip("\n").split("\t", 1)
            id_map[sanitized] = original

    with open(csv_path, newline="") as fh, open(out_path, "w", newline="") as out:
        reader = csv.reader(fh)
        writer = csv.writer(out)
        writer.writerow(next(reader))  # header

        for row in reader:
            query = row[0]
            if query.strip():
                sanitized = query.replace("|", "_")
                if sanitized not in id_map:
                    sys.exit(
                        f"ERROR: sanitized id '{sanitized}' (from QUERY "
                        f"'{query}') has no entry in the id map -- the query "
                        "fasta used to build the map must not be the exact "
                        "fasta TreeGrafter ran on."
                    )
                row[0] = id_map[sanitized]
            writer.writerow(row)


if __name__ == "__main__":
    main()
