


process ORTHOLOGER {

    tag "Logging $project_id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/orthologer' }"
    

    input:
    tuple val(meta), path (odbdata)
    tuple val(meta), path (odbwork)
    val path_to_work
    val project_id 

    output:
    tuple val(meta), path ("$odbwork/Cluster/odbproj.og"), emit: loger
    path ("$odbwork/versions.yml")               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    cp -r $odbdata/* /odbdata
    cd /odbwork
    setup_odb.sh
    orthologer -c import
    cd Sequences/
    capitalize_files.py
    cd ../
    orthologer \\
        -c run \\
        ${args}
    cp -r Cluster/ /tmp/*/odbwork/ 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthologer: \$(./orthologer.sh -v)
    END_VERSIONS
    cp versions.yml /tmp/*/odbwork/
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthologer: \$(./orthologer.sh -v)
    END_VERSIONS
    """

}
