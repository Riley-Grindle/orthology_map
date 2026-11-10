
process POST_PROC {
    tag "Processing orthology tool outputs..."
    //label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/post_proc' }"

    input:
    tuple val(meta), path (outs)
    tuple val(meta), path (fastas)
    val fasta_file
    path ensembl_data    // absolute path or null — not staged, accessed directly on shared fs
    val taxa_db         // absolute path or null — not staged, accessed directly on shared fs
    path tree_grafter_raw       // TREEGRAFTER's own un-split "*.out.txt" (ch_tree) — query -> (family, subfamily), used only to compute NEAREST_REF_NODE_DISTANCE below
    path panther_db_data        // same panther_db_data dir TREEGRAFTER was given (params.sup_data_panther) — per-family Newick trees + leaf annotations
    val ref_panther_prefixes    // same comma-separated PANTHER gene-ID prefix list PANTHER_API uses, e.g. "MGI"

    output:
    tuple val(meta), path ("*ortho_f_std.tsv"), emit: ortho_f
    tuple val(meta), path ("*ortho_l_std.tsv"), emit: ortho_l
    path ("*eggnog_std.tsv"),  emit: eggnog
    path ("*tree_std.tsv"),    emit: tree

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    FILE="\$(basename $fasta_file)"
    PREFIX="\$(echo "\$FILE" | rev | cut -d. -f2- | rev)"
    
    if ! [[ -e "./outs/ortho_f/blank.txt" ]]; then
        Rscript ${projectDir}/bin/ortho_f_post.R ./outs/ortho_f/*/Orthologues/Orthologues_\${PREFIX}
    fi
        
    # ORTHOLOGER sanitizes the query label (replaces any char outside
    # [A-Za-z0-9_] with '_') before writing it into the .og cluster file,
    # so this must match that sanitization or ortho_l_post.R finds 0 rows.
    # (tr -c also translates the newline `cut` terminates its output with,
    # so that newline must be stripped -- via the \$(...) capture -- before
    # piping into tr, or it becomes a spurious trailing '_'.)
    ORTHO_L_RAW="\$(echo "\$FILE" | cut -d. -f1)"
    ORTHO_L_PREFIX="\$(printf '%s' "\$ORTHO_L_RAW" | tr -c 'A-Za-z0-9_' '_')"
    Rscript ${projectDir}/bin/ortho_l_post.R ./outs/ortho_l/ "\${ORTHO_L_PREFIX}"
    Rscript ${projectDir}/bin/eggnog_post.R ./outs/egg
    Rscript ${projectDir}/bin/tree_post.R ./outs/tree
    mkdir headers
    for file in ./odbdata/*; do
        FILE_NAME="\$(basename "\$file")"; filename_no_extension="\${FILE_NAME%.*}_headers.txt"
        grep ">" "\$file" | grep -o '^[^[:space:]]*' > "headers/"\$filename_no_extension""
    done
    Rscript ${projectDir}/bin/hexadecimal_correction.R ./odbdata ./headers/ ./
    
    mv ./tree_std.csv input_tree.csv
    build_treegrafter_id_map.py "$fasta_file" tree_id_map.tsv
    fix_tree_query_ids.py input_tree.csv tree_id_map.tsv input_tree.fixed.csv
    mv input_tree.fixed.csv input_tree.csv
    if [ -n "$ensembl_data" ]; then
        ensembl_id_2_gene_symbl.py ./input_tree.csv /sup_data/prefix_2_file.json $ensembl_data /sup_data/prefix_2_species.json
        paste -d"," input_tree.csv gene_symbols_ensembl.txt species_names.txt > tree_std.csv
    else
        cp input_tree.csv tree_std.csv
    fi
    if [ -n "$taxa_db" ]; then
        search_taxa.sh tree_std.csv $taxa_db
        echo "SPECIES" >> tmp.tsv
        cat taxa.tsv >> tmp.tsv
        sed 's/^/"/' tmp.tsv | sed 's/\$/"/' > taxa_col.tsv
        cut -d, -f 1-3 tree_std.csv > editing_tree.csv
        paste -d"," editing_tree.csv taxa_col.tsv > tree_std.csv
    fi
    # Appends NEAREST_REF_NODE_DISTANCE/NEAREST_REF_GENE_ID/NEAREST_REF_GENE_NAME:
    # for every "NO MATCH" row, how many PANTHER12 gene-tree edges separate
    # that query's TreeGrafter subfamily from the nearest reference-species
    # gene, and which gene that actually is (a distance with no id attached
    # to it isn't reviewable -- NA for real matches and continuation rows).
    # nearest_ref_node_distance.py writes its own already-quoted CSV
    # fragment (same QUOTE_ALL convention as tree_std.csv), so this just
    # pastes it on directly.
    if [ -n "$panther_db_data" ]; then
        nearest_ref_node_distance.py tree_std.csv tree_id_map.tsv "$tree_grafter_raw" "$panther_db_data" "${ref_panther_prefixes}" distance_cols.csv
        paste -d"," tree_std.csv distance_cols.csv > tree_std_with_distance.csv
        mv tree_std_with_distance.csv tree_std.csv
    fi
    cut -d"," -f2 eggnog_std.csv | sed "s/^.//" | sed "s/.\$//" > match_ids.txt
    mv eggnog_std.csv input_eggnog.csv
    string_search.py /sup_data/string_db.json ./match_ids.txt
    sed 's/^/"/' gene_symbols.txt | sed 's/\$/"/' > gene_symbols_01.txt
    sed 's/^/"/' taxa_ids.txt | sed 's/\$/"/' > taxa_ids_01.txt
    paste -d "," input_eggnog.csv gene_symbols_01.txt taxa_ids_01.txt > eggnog_std.csv
    
    for file in *std.csv; do
        sed 's/^.//' \${file} | sed 's/","/\t/g' | sed 's/.\$//' > ${meta.id}.\${file::-4}.tsv
    done
    # Harmless cache-busting comment, v2: nearest_ref_node_distance.py now
    # also reports which reference gene a NO-MATCH row's distance is
    # measured to (NEAREST_REF_GENE_ID/NAME), not just the hop count --
    # bin/ script changes alone don't invalidate -resume's cache, only a
    # change to this process's own script body does.
    """

}
