

process ORTHOFINDER {
    
    tag "Orthofinding $project_id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/orthofinder' }"

    input:
    path fasta
    val project_id

    output:
    path ("$fasta/OrthoFinder/Results_*"), emit: ortho_f
    path ("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    orthofinder \\
        -f $fasta \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(echo \$(orthofinder version) | sed "s/orthofinder, version //g" )
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """
}

