

process COMBINE_JSON {
    
    tag "${project_id}:${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/orthofinder' }"

    input:
    tuple val(meta), path (fasta)
    val project_id

    output:
    tuple val(meta), path ("$fasta/OrthoFinder/Results_*"), emit: ortho_f
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
        orthofinder: \$(echo \$(orthofinder --version) | sed "s/orthofinder, version //g" )
    END_VERSIONS
    """
}

