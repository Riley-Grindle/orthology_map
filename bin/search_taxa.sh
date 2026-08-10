#!/usr/bin/bash
#
# For each data row of $1's 4th (SPECIES) column, look up a taxon in $2 and
# append exactly one line to taxa.tsv. Two independent alignment bugs fixed
# here, both of which corrupt the paste that reassembles tree_std.csv
# downstream in post_proc.nf:
#
# 1. The original version only appended when `grep | head -n1 | cut -f1`
#    produced output, so a query that found no match in $2 wrote nothing at
#    all -- silently desyncing taxa.tsv's line count from the input.
#    Always writing a line (falling back to "NA") keeps every row aligned
#    regardless of whether the lookup succeeds.
# 2. `cut -d, -f 4 "$1"` included $1's own header row, extracting the
#    literal string "SPECIES" as if it were a real value to look up.
#    post_proc.nf already prepends its own "SPECIES" header
#    (`echo "SPECIES" >> tmp.tsv`) before appending this script's output, so
#    processing $1's header here produced one extra leading line and shifted
#    every real data row's taxon by one position. tail -n +2 skips it.

tail -n +2 "$1" | cut -d, -f 4 | sed 's/"//g' | while IFS= read -r line; do
    result=""
    if [ "$line" != "NA" ]; then
        result="$(grep "$line" "$2" | head -n 1 | cut -f 1)"
    fi
    if [ -n "$result" ]; then
        echo "$result" >> taxa.tsv
    else
        echo "NA" >> taxa.tsv
    fi
done
