
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
    FILE=$fasta_file
    PREFIX="\$(echo "\$FILE" | cut -d. -f1)"
    Rscript /rscripts/ortho_f_post.R ./outs/ortho_f/*/Orthologues/Orthologues_\${PREFIX}
    Rscript /rscripts/ortho_l_post.R ./outs/ortho_l/ "\${PREFIX}"
    Rscript /rscripts/eggnog_post.R ./outs/egg
    Rscript /rscripts/tree_post.R ./outs/tree
    mkdir /headers
    for file in ./odbdata/*; do
        FILE_NAME="\$(basename "\$file")"; filename_no_extension="\${FILE_NAME%.*}_headers.txt"
        grep ">" "\$file" | grep -o '^[^[:space:]]*' > "/headers/"\$filename_no_extension""
    done
    Rscript /rscripts/hexadecimal_correction.R ./odbdata /headers/ ./
    mkdir std_outs
    mv *std.csv std_outs/
    """

}
