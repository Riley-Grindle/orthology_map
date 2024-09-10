


process DIAMOND {
    tag "Making Diamond db for $project_id using: $fasta"
    label "process_high"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/diamond' }"

    input:
    tuple val(meta), path (fasta)


    output:
    path ("*.dmnd"), emit: egg
    //path ("versions.yml")            , emit: versions

    script:
    def args = task.ext.args  ?: ''
    """
    diamond makedb --in $fasta -d $meta

    """

    stub:
    """
    touch test.txt
    """
}

