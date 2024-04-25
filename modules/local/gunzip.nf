
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
    gtf_name=$gtf
    fasta_name=$fasta
    ext="\${fasta_name##*.}"
   
    if [[ "\$gtf_name" == *".gz"* ]]; then
        gunzip $gtf
        cp *.gtf ${meta.id}.gtf
    else
        cp $gtf ${meta.id}.gtf
    fi

    if [[ "\$fasta_name" == *".gz"* ]]; then
        gunzip $fasta
        cp *.\$ext ${meta.id}.fa
    else
        cp *.\$ext ${meta.id}.fa
    fi
    """

    stub:
    """
    touch $gtf
    touch $fasta
    """


}
