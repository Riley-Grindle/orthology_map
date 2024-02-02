process COLLECT_CHUNKS {
    tag "Combining Outputs"
    label "process_low"
 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/python:3.9--1' }"
 
    input:
    path ortho_text

    output:
    path "orthologs_final.txt"

    script:
    """
    touch orthologs_final.txt
    cat $ortho_text >> orthologs_final.txt
    """

}
