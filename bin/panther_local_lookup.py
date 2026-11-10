#!/usr/bin/env python3
"""
Reconstructed 2026-11-10 after the original bin/panther_local_lookup.py was
lost in an infrastructure incident (JuiceFS filesystem destroyed and
restored from a metadata dump that predates this file's creation). Rebuilt
from the calling contract documented in bin/post_panther.sh's comments and
reverse-engineered from bin/tree_post.R's parser (both of which survived),
plus the real PANTHER12.0 "sequence classification" bulk files fetched
fresh from data.pantherdb.org. Verified against the one concrete example
left behind in post_panther.sh's comments: PTHR23333:SF24 -> MGI=3042273.

Usage:
    panther_local_lookup.py SUBFAMILIES_FILE CLASSIFICATION_DIR REF_PREFIXES OUT_DIR

SUBFAMILIES_FILE: one PANTHER id per line, in the same order as
    post_panther.sh's sequence_headers.txt (1:1 positional correspondence).
    Each line is either "PTHR#####:SF##" (subfamily resolved) or just
    "PTHR#####" (TreeGrafter only resolved to family level).
CLASSIFICATION_DIR: directory containing one file per PANTHER gene-ID
    prefix (e.g. "MGI", "ZFIN", "HGNC" -- exactly the prefix strings in
    assets/panther_species_prefixes.json), each a tab-delimited PANTHER12.0
    sequence-classification file (column 1 = gene identifier
    "SPECIES|PREFIX=id|UniProtKB=...", column 3 = "PTHR#####[:SF##]").
REF_PREFIXES: comma-separated prefixes to search (e.g. "MGI" or
    "MGI,ZFIN"); empty means no configured reference has a known prefix,
    so every output file is empty (matches pre-existing fallback behavior).
OUT_DIR: directory to write genes1.txt, genes2.txt, ... into, one per line
    of SUBFAMILIES_FILE, each line "gene_id=<value>" (empty file = no
    match). This exact "gene_id=" (8-char) prefix is required by
    tree_post.R, which locates match lines by checking
    substr(line, 1, 7) == "gene_id" and reads the value from substr(line, 9, ...).
"""

import os
import sys


def load_index(path):
    """subfamily/family id (column 3, exact) -> list of normalized gene ids (column 1, normalized)."""
    index = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            gene_field, pthr_id = fields[0], fields[2]
            if not pthr_id.startswith("PTHR"):
                continue
            index.setdefault(pthr_id, []).append(gene_field)
    return index


def normalize_gene_id(gene_field, prefix):
    """
    gene_field looks like "MOUSE|MGI=MGI=3042273|UniProtKB=Q9CZ44" (species |
    PREFIX=native_id | UniProtKB=accession, pipe-delimited, only the middle
    field is a prefix match). Some PANTHER12 rows use a different DB label
    in the middle field (Ensembl, Gene, GeneID, ...) instead of the
    species' PANTHER-tree prefix -- those genes just aren't reachable via
    this prefix-keyed lookup, same recall gap as before this rebuild.

    The native id for MGI/RGD-style organisms is itself already
    "<PREFIX>:<num>" (": " swapped for "=" by PANTHER), producing a doubled
    "MGI=MGI=3042273" -- collapse exactly one such duplicate layer so the
    emitted value matches the single "MGI=3042273" form confirmed against
    production (post_panther.sh's own PTHR23333:SF24 example). Organisms
    whose native id does NOT repeat the prefix (ZFIN, FlyBase, WormBase,
    SGD, CGD, PomBase, dictyBase, TAIR, HGNC, Xenbase) pass through with
    exactly one "<PREFIX>=" layer, unchanged.
    """
    parts = gene_field.split("|")
    if len(parts) < 2:
        return None
    mid = parts[1]
    marker = prefix + "="
    if not mid.startswith(marker):
        return None
    rest = mid[len(marker):]
    if rest.startswith(marker):
        rest = rest[len(marker):]
    return f"{prefix}={rest}"


def lookup(query_id, index):
    """Exact subfamily match; family-only ids also pull in every subfamily under that family."""
    matches = list(index.get(query_id, []))
    if ":" not in query_id:
        family_prefix = query_id + ":"
        for key, values in index.items():
            if key.startswith(family_prefix):
                matches.extend(values)
    return matches


def main():
    if len(sys.argv) != 5:
        sys.stderr.write(
            "usage: panther_local_lookup.py SUBFAMILIES_FILE CLASSIFICATION_DIR REF_PREFIXES OUT_DIR\n"
        )
        sys.exit(1)

    subfamilies_file, classification_dir, ref_prefixes, out_dir = sys.argv[1:5]

    prefixes = [p.strip() for p in ref_prefixes.split(",") if p.strip()]

    indexes = {}
    for prefix in prefixes:
        path = os.path.join(classification_dir, prefix)
        if not os.path.isfile(path):
            sys.stderr.write(
                f"WARNING: no classification file for prefix '{prefix}' at {path}; "
                "treating as no matches for this prefix\n"
            )
            continue
        indexes[prefix] = load_index(path)

    os.makedirs(out_dir, exist_ok=True)

    with open(subfamilies_file, "r", encoding="utf-8", errors="replace") as f:
        query_ids = [line.rstrip("\n") for line in f]

    for i, query_id in enumerate(query_ids, start=1):
        out_path = os.path.join(out_dir, f"genes{i}.txt")
        lines = []
        for prefix in prefixes:
            index = indexes.get(prefix)
            if index is None:
                continue
            for gene_field in lookup(query_id, index):
                value = normalize_gene_id(gene_field, prefix)
                if value is not None:
                    lines.append(f"gene_id={value}")
        with open(out_path, "w", encoding="utf-8") as out:
            if lines:
                out.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
