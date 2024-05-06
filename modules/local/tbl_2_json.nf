

process TBL_2_JSON {
    
    tag "${meta.tool}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/tool_vote:latest' }"

    input:
    tuple val(meta), path (outfile)

    output:
    tuple val(meta), path ("*.query_2_matches.json"), emit: json
    //path ("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    tbl_2_json.py \\
        $outfile \\
        ${args} \\
        ${meta.tool}

    """

    stub:
    """
    touch $outfile
    """
}

