#!/usr/bin/env python3
"""
For every row in tree_std.csv whose MATCH is "NO MATCH", compute how many
PANTHER12 gene-tree edges separate that query's TreeGrafter subfamily from
the nearest leaf belonging to one of this run's configured reference
species, and append the result as three new columns:
  - NEAREST_REF_NODE_DISTANCE: the hop count itself
  - NEAREST_REF_GENE_ID: which reference gene that distance is measured
    to, in the same "<PREFIX>:<value>" format the MATCH column itself uses
    for a real hit (e.g. "MGI:MGI=107792"), so the two columns are directly
    comparable/joinable
  - NEAREST_REF_GENE_NAME: that gene's PAINT-curated protein name (e.g.
    "LIM/homeobox protein Lhx5"), for a human reading the file without a
    separate id lookup

A distance with no accompanying id is useless for review -- you can't tell
whether "11" means "practically the same gene" or "a barely-related
paralog" without knowing what's on the other end of it.

Rows that aren't "NO MATCH" (a real gene_id match, or a blank continuation
row -- see tree_post.R) get "NA" in all three columns: a real match already
has a reference ortholog, so a tree distance adds nothing.

Data sources, all genuine PANTHER12 data TreeGrafter itself ships with (no
live API calls):
  - tree_grafter_raw: the TREEGRAFTER process's own "*.out.txt", mapping
    TreeGrafter's sanitized query id -> (family, subfamily)
  - panther_db_data: PANTHER12_Tree_MSF/<FAMILY>.newick (per-family gene
    tree) and PANTHER12_PAINT_Annotations/PANTHER12_leaf_GO_IBDannotations.tab
    (leaf AN id -> subfamily + species|gene, used both to find which leaves
    are reference-species leaves and which leaves belong to a subfamily)
  - tree_id_map.tsv: already built earlier in post_proc.nf by
    build_treegrafter_id_map.py (original fasta header -> TreeGrafter's
    sanitized id) -- reused here (reversed) instead of recomputing the
    sanitization regex, so there is one source of truth for it.

Method: parse each needed family's Newick tree into a plain node graph
(branch lengths ignored -- this is topological/node distance, not branch
length), identify which leaves belong to a configured reference species and
which leaves belong to the query's subfamily, then multi-source BFS from
all reference-species leaves at once. The reported distance is the min hop
count over all of the subfamily's member leaves.
"""

import re
import sys
from pathlib import Path
from collections import deque, defaultdict
import csv

# PANTHER's gene-id prefix -> the species code used as the first "|"
# delimited token in panther_db_data's bulk/annotation files (confirmed
# against the PANTHER12.0 per-organism classification files already fetched
# for panther_local_lookup.py: PTHR12.0_mouse starts every row "MOUSE|...",
# _zebrafish "DANRE|...", _human "HUMAN|...", _rat "RAT|..."). Only prefixes
# actually verified are listed; an unrecognized prefix is skipped (NA
# output) with a warning rather than treated as a hard failure, since this
# column is a diagnostic enrichment, not something the rest of the pipeline
# depends on.
PREFIX_TO_SPECIES_CODE = {
    "MGI": "MOUSE",
    "ZFIN": "DANRE",
    "HGNC": "HUMAN",
    "RGD": "RAT",
}
SPECIES_CODE_TO_PREFIX = {v: k for k, v in PREFIX_TO_SPECIES_CODE.items()}

SF_RE = re.compile(r"SF:(PTHR\d+(?::SF\d+)?)")
NAME_RE = re.compile(r"NAME:([^;]*)")


def normalize_gene_id(gene_info_raw, prefix):
    # gene_info_raw looks like "MOUSE|MGI=MGI=107792|UniProtKB=P61375" --
    # find the "<prefix>=<value>" segment and replace only the first "="
    # with ":" to match the MATCH column's own "<prefix>:<value>"
    # convention (same transform as panther_local_lookup.py uses when
    # building that column in the first place).
    for part in gene_info_raw.split("|"):
        if part.startswith(prefix + "="):
            key, _, value = part.partition("=")
            return f"{key}:{value}"
    return None


def load_sanitized_lookup(tree_id_map_path):
    """original fasta header -> TreeGrafter's sanitized id (reverse of the
    sanitized -> original map build_treegrafter_id_map.py already built)."""
    original_to_sanitized = {}
    with open(tree_id_map_path) as fh:
        for line in fh:
            sanitized, original = line.rstrip("\n").split("\t", 1)
            original_to_sanitized[original] = sanitized
    return original_to_sanitized


def load_query_to_subfamily(tree_grafter_raw_path):
    """sanitized query id -> (family, subfamily_or_None)"""
    mapping = {}
    with open(tree_grafter_raw_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            qid, family, annot = parts[0], parts[1], parts[2]
            m = SF_RE.search(annot)
            subfamily = m.group(1) if (m and ":SF" in m.group(1)) else None
            mapping[qid] = (family, subfamily)
    return mapping


class NewickParser:
    """Minimal Newick parser -> adjacency dict {node_name: set(neighbors)}.
    PANTHER12 format: every node (leaf and internal) is named, branch
    lengths follow ':', no NHX comments."""

    def __init__(self, text):
        self.text = text
        self.pos = 0
        self.edges = defaultdict(set)
        self.leaves = set()
        self._anon = 0

    def parse(self):
        self._parse_clade()
        return self.edges, self.leaves

    def _name_and_len(self):
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in ",();:":
            self.pos += 1
        name = self.text[start:self.pos]
        if self.pos < len(self.text) and self.text[self.pos] == ":":
            self.pos += 1
            while self.pos < len(self.text) and self.text[self.pos] not in ",();":
                self.pos += 1
        if not name:
            self._anon += 1
            name = f"__anon{self._anon}"
        return name

    def _parse_clade(self):
        if self.text[self.pos] == "(":
            self.pos += 1
            children = []
            while True:
                children.append(self._parse_clade())
                if self.text[self.pos] == ",":
                    self.pos += 1
                    continue
                self.pos += 1  # ')'
                break
            me = self._name_and_len()
            for c in children:
                self.edges[me].add(c)
                self.edges[c].add(me)
            return me
        name = self._name_and_len()
        self.leaves.add(name)
        return name


def bfs_multi_source(edges, sources):
    """Returns (dist, origin): origin[node] is which source leaf the
    shortest path to `node` actually came from, so a per-leaf distance can
    be traced back to a specific reference gene, not just a number."""
    dist = {s: 0 for s in sources}
    origin = {s: s for s in sources}
    q = deque(sources)
    while q:
        node = q.popleft()
        for nxt in edges.get(node, ()):
            if nxt not in dist:
                dist[nxt] = dist[node] + 1
                origin[nxt] = origin[node]
                q.append(nxt)
    return dist, origin


def main():
    if len(sys.argv) != 7:
        sys.exit(
            f"Usage: {sys.argv[0]} <tree_std.csv> <tree_id_map.tsv> "
            f"<tree_grafter_raw.out.txt> <panther_db_data_dir> "
            f"<ref_panther_prefixes_csv> <output.tsv>"
        )

    tree_std_csv, tree_id_map_path, raw_path, db_dir, ref_prefixes_arg, out_path = sys.argv[1:7]
    db_dir = Path(db_dir)
    newick_dir = db_dir / "PANTHER12_Tree_MSF"
    leaf_annot_path = db_dir / "PANTHER12_PAINT_Annotations" / "PANTHER12_leaf_GO_IBDannotations.tab"

    ref_species_codes = set()
    for prefix in (p for p in ref_prefixes_arg.split(",") if p):
        code = PREFIX_TO_SPECIES_CODE.get(prefix)
        if code is None:
            print(
                f"WARNING: no known PANTHER species code for prefix '{prefix}' -- "
                f"NEAREST_REF_NODE_DISTANCE will be NA wherever only this "
                f"prefix would apply. Add it to PREFIX_TO_SPECIES_CODE in "
                f"{sys.argv[0]} if needed.",
                file=sys.stderr,
            )
            continue
        ref_species_codes.add(code)

    original_to_sanitized = load_sanitized_lookup(tree_id_map_path)
    query_to_subfamily = load_query_to_subfamily(raw_path)

    # Read tree_std.csv rows up front so we know exactly which families we
    # need before paying for the ~1M-line leaf annotation file.
    rows = []
    with open(tree_std_csv, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for row in reader:
            rows.append(row)

    needed = []  # (row_index, family, subfamily)
    for i, row in enumerate(rows):
        if len(row) < 2:
            needed.append((i, None, None))
            continue
        query, match = row[0].strip(), row[1].strip()
        if not query or match != "NO MATCH":
            needed.append((i, None, None))
            continue
        sanitized = original_to_sanitized.get(query)
        info = query_to_subfamily.get(sanitized) if sanitized else None
        if info is None or info[1] is None:
            needed.append((i, None, None))
            continue
        family, subfamily = info
        needed.append((i, family, subfamily))

    families_needed = {f for _, f, sf in needed if f is not None}

    family_leaf_info = {f: {} for f in families_needed}
    if families_needed and ref_species_codes and leaf_annot_path.is_file():
        with open(leaf_annot_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                an_full, annot, gene_info = parts[0], parts[1], parts[2]
                family, _, an_id = an_full.partition(":")
                if family not in family_leaf_info:
                    continue
                m = SF_RE.search(annot)
                subfamily = m.group(1) if (m and ":SF" in m.group(1)) else None
                name_m = NAME_RE.search(annot)
                name = name_m.group(1).strip() if name_m else None
                species_code = gene_info.split("|", 1)[0]
                family_leaf_info[family][an_id] = (species_code, subfamily, gene_info, name)

    dist_per_leaf_cache = {}  # family -> (dist_map, origin_map) or None if unusable
    sf_nearest_cache = {}  # (family, subfamily) -> (distance, nearest_an) or (None, None)

    def get_family_dist_map(family):
        if family in dist_per_leaf_cache:
            return dist_per_leaf_cache[family]
        leaf_info = family_leaf_info.get(family, {})
        ref_leaves = {an for an, info in leaf_info.items() if info[0] in ref_species_codes}
        if not ref_leaves:
            dist_per_leaf_cache[family] = None
            return None
        newick_path = newick_dir / f"{family}.newick"
        try:
            with open(newick_path) as fh:
                text = fh.read().strip()
        except FileNotFoundError:
            dist_per_leaf_cache[family] = None
            return None
        edges, leaves = NewickParser(text).parse()
        valid_ref_leaves = ref_leaves & leaves
        if not valid_ref_leaves:
            dist_per_leaf_cache[family] = None
            return None
        result = bfs_multi_source(edges, valid_ref_leaves)
        dist_per_leaf_cache[family] = result
        return result

    def get_subfamily_nearest(family, subfamily):
        """-> (distance, gene_id_or_None, gene_name_or_None)"""
        key = (family, subfamily)
        if key in sf_nearest_cache:
            return sf_nearest_cache[key]
        result = get_family_dist_map(family)
        if result is None:
            sf_nearest_cache[key] = (None, None, None)
            return sf_nearest_cache[key]
        dist_map, origin_map = result
        leaf_info = family_leaf_info.get(family, {})
        members = [an for an, info in leaf_info.items() if info[1] == subfamily]
        candidates = [(dist_map[an], origin_map[an]) for an in members if an in dist_map]
        if not candidates:
            sf_nearest_cache[key] = (None, None, None)
            return sf_nearest_cache[key]
        distance, nearest_an = min(candidates, key=lambda c: c[0])
        nearest_species, _, nearest_gene_info, nearest_name = leaf_info[nearest_an]
        prefix = SPECIES_CODE_TO_PREFIX.get(nearest_species)
        gene_id = normalize_gene_id(nearest_gene_info, prefix) if prefix else None
        sf_nearest_cache[key] = (distance, gene_id, nearest_name)
        return sf_nearest_cache[key]

    results = [(None, None, None)] * len(rows)
    for i, family, subfamily in needed:
        if family is None:
            continue
        results[i] = get_subfamily_nearest(family, subfamily)

    with open(out_path, "w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n", quoting=csv.QUOTE_ALL)
        writer.writerow(["NEAREST_REF_NODE_DISTANCE", "NEAREST_REF_GENE_ID", "NEAREST_REF_GENE_NAME"])
        for distance, gene_id, name in results:
            writer.writerow([
                distance if distance is not None else "NA",
                gene_id if gene_id is not None else "NA",
                name if name is not None else "NA",
            ])


if __name__ == "__main__":
    main()
