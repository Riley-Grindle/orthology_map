#!/bin/bash

# $1: comma-separated PANTHER gene-ID prefixes for the reference species in
# this run (e.g. "MGI" or "MGI,ZFIN"), resolved upstream in orthologymap.nf
# from each reference's NCBI taxon ID via assets/panther_species_prefixes.json.
# Empty means no known prefix for the configured reference(s): every query
# below then matches nothing (genes/geneN.txt all empty), same as before.
REF_PREFIXES="$1"

# PANTHER12.0 per-organism "sequence classification" bulk files, pre-fetched
# from data.pantherdb.org (the official archive; not the live treeinfo API)
# onto JuiceFS so every Float task can read them directly with no network
# call at run time. See bin/panther_local_lookup.py for why this replaced
# querying the live API (version skew + reliability).
CLASSIFICATION_DIR="/mnt/jfs/nextflow/panther12_classifications"

# TreeGrafter's raw ".out" rows are tab-delimited:
# <query_id>\t<family>\t<annotation>, where <annotation> looks like
# "SF:<family>:<subfamily>;NAME:<curated name>;GO-Annotations:...;" (or just
# "SF:<family>;..." when TreeGrafter could only resolve to the family
# level). Filtering/extracting by tab-delimited field (awk -F'\t'), not a
# single free-floating regex over the whole line: the previous
# `grep -o "^.*.PTHR.*.SF[^;]*"` could re-anchor on a literal "SF"
# substring embedded in the curated NAME itself (e.g. "NSFL1 cofactor p47"
# contains "SF"), silently extending the "match" past the real subfamily id
# and gluing trailing ";NAME:..." text onto it. That corrupted id then
# never matched the clean key in panther_local_lookup.py's index, producing
# a false "NO MATCH" for a gene that DOES have a real reference-species
# ortholog in that exact subfamily (confirmed: PTHR23333:SF24 has a mouse
# member, MGI=3042273, in the PANTHER12 bulk data, but a zebrafish query
# grafted onto it came back with zero genes purely because of this bug).
awk -F'\t' '$2 ~ /^PTHR[0-9]+$/ && $3 ~ /^SF:/' *.out > cleaned.txt

sed 's/\(PTHR.*\)//' cleaned.txt > sequence_headers.txt

touch "orthologs.txt"

awk -F'\t' '{ann = $3; sub(/^SF:/, "", ann); sub(/;.*$/, "", ann); print ann}' cleaned.txt > subfamilies.txt

mkdir "./genes"

PY=python3
command -v "$PY" >/dev/null 2>&1 || PY=python
command -v "$PY" >/dev/null 2>&1 || {
    echo "ERROR: neither python3 nor python found on PATH" >&2
    exit 1
}
LOOKUP_SCRIPT="$(command -v panther_local_lookup.py)"
if [ -z "$LOOKUP_SCRIPT" ]; then
    echo "ERROR: panther_local_lookup.py not found on PATH" >&2
    exit 1
fi

"$PY" "$LOOKUP_SCRIPT" subfamilies.txt "$CLASSIFICATION_DIR" "$REF_PREFIXES" "./genes"

count=1
while IFS= read -r id; do

    gene_file="./genes/genes$count.txt"

    echo "$id" >> "orthologs.txt"
    echo -e "_________________________\n" >> "orthologs.txt"
    cat "$gene_file" >> "orthologs.txt"

    echo -e "\n_______________________________________________________________________\n" >> "orthologs.txt"

    count=$((count+1))

done < sequence_headers.txt
