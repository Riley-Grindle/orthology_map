
process GUNZIP {
    tag "Splitting inputs into channels"
    label "process_low"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(gtf), path(fasta)

    output:
    tuple val(meta), path("*.gtf"), emit: gtf
    tuple val(meta), path("*.fa"), emit: fasta

    script:
    """
    cp $gtf ${meta.id}.gtf
    cp $fasta ${meta.id}.fa
    """

    stub:
    """
    touch $gtf
    touch $fasta
    """


}
