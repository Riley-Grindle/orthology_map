#!/usr/bin/env python3
"""
Add a per-match support tally and a designated best_match to each query
gene in voted_orthologs.json, without dropping any evidence.

Up to this point "voting" has really been a union: every match any tool
proposed for a query gene is kept, tagged with which (species, tool) pairs
back it. This script is the first place an actual inclusion/confidence
decision happens -- everything stays, but each match now records how many
independent tools support it, matches are ranked, and the top-ranked match
per gene is exposed as best_match for downstream table-building code to use.

Run after species_layer.py, which reshapes vote.py's (species, tool) tuples
into {match: {species: [tools]}}.
"""

import sys
import json
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("json_path")
    parser.add_argument(
        "--min-tool-support",
        type=int,
        default=2,
        help="matches backed by fewer distinct tools than this are tagged low-confidence",
    )
    args = parser.parse_args()

    with open(args.json_path) as fh:
        aggregate = json.load(fh)

    for query, matches in aggregate.items():
        ranked = []
        for symbol, species_map in matches.items():
            tools = set()
            for tool_list in species_map.values():
                tools.update(tool_list)
            ranked.append(
                {
                    "symbol": symbol,
                    "species": species_map,
                    "tool_count": len(tools),
                    "species_count": len(species_map),
                }
            )

        ranked.sort(key=lambda m: (-m["tool_count"], -m["species_count"], m["symbol"]))

        new_matches = {}
        for i, m in enumerate(ranked, start=1):
            new_matches[m["symbol"]] = {
                "species": m["species"],
                "tool_count": m["tool_count"],
                "species_count": m["species_count"],
                "rank": i,
                "confidence": "high" if m["tool_count"] >= args.min_tool_support else "low",
            }

        aggregate[query] = {
            "matches": new_matches,
            "best_match": ranked[0]["symbol"] if ranked else None,
        }

    with open(args.json_path, "w") as fh:
        json.dump(aggregate, fh, indent=4)


if __name__ == "__main__":
    main()
