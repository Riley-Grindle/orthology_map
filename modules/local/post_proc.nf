
process POST_PROC {
    tag "Processing orthology tool outputs..."
    label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/post_proc' }"

    input:
    path outs
    path fastas
    val fasta_file

    output:
    path ("./std_outs/")

    script:
    """
    FILE="\$(basename $fasta_file)"
    PREFIX="\$(echo "\$FILE" | rev | cut -d. -f2- | rev)"
    if ! [[ -e "./outs/ortho_f/blank.txt" ]]; then
        Rscript /rscripts/ortho_f_post.R ./outs/ortho_f/*/Orthologues/Orthologues_\${PREFIX}
    fi
    ORTHO_L_PREFIX="\$(echo "\$FILE" | cut -d. -f1)"
    Rscript /rscripts/ortho_l_post.R ./outs/ortho_l/ "\${ORTHO_L_PREFIX}"
    Rscript /rscripts/eggnog_post.R ./outs/egg
    Rscript /rscripts/tree_post.R ./outs/tree
    mkdir /headers
    for file in ./odbdata/*; do
        FILE_NAME="\$(basename "\$file")"; filename_no_extension="\${FILE_NAME%.*}_headers.txt"
        grep ">" "\$file" | grep -o '^[^[:space:]]*' > "/headers/"\$filename_no_extension""
    done
    Rscript /rscripts/hexadecimal_correction.R ./odbdata /headers/ ./
    Rscript /rscripts/biomart.R tree_std.csv
    cut -d"," -f2 eggnog_std.csv | sed "s/^.//" | sed "s/.\$//" > match_ids.txt
    string_search.py /sup_data/string_db.json ./match_ids.txt
    sed 's/^/"/' gene_symbols.txt | sed 's/\$/"/' > gene_symbols_01.txt
    paste -d "," eggnog_std.csv gene_symbols_01.txt
    mkdir std_outs
    mv *std.csv std_outs/
    """

}
