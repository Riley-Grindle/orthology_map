#!/usr/bin/env python3
"""
Build a lookup from TreeGrafter's sanitized query ids back to the original
FASTA header ids.

treeGrafter.pl sanitizes every non-word character in the first
whitespace-delimited token of each header to "_" before using it as the
query id in its own output (`$queryid =~ s/[^\\w]/\\_/g;`). Our query ids
are TransDecoder's "GENE.<gene>~~<transcript>.<pN>" headers, so "~~" and "."
both collapse to "_" -- and because "_" is itself a word character, any
underscore already present in the gene or transcript id survives untouched
and becomes indistinguishable from an injected separator once sanitized.

Pattern-matching the sanitized string back into gene/transcript components
(the previous approach, in tree_post.R) is therefore ambiguous whenever an
id naturally contains an underscore. This script instead replays the exact
same forward transformation on the true headers, so the post-processing
step can do an exact dictionary lookup instead of guessing.
"""

import re
import sys

SANITIZE = re.compile(r"[^A-Za-z0-9_]")


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <query.fasta> <output.tsv>")

    fasta_path, out_path = sys.argv[1], sys.argv[2]

    mapping = {}
    with open(fasta_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if not line.startswith(">"):
                continue
            original = line[1:].split()[0]
            sanitized = SANITIZE.sub("_", original)
            if sanitized in mapping and mapping[sanitized] != original:
                sys.exit(
                    f"ERROR: '{mapping[sanitized]}' and '{original}' both "
                    f"sanitize to '{sanitized}' -- TreeGrafter cannot tell "
                    "them apart, so the id map would be ambiguous for this pair."
                )
            mapping[sanitized] = original
            out.write(f"{sanitized}\t{original}\n")


if __name__ == "__main__":
    main()
