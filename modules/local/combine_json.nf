

process COMBINE_JSON {
    
    tag "${meta.tool}"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/tool_vote:latest' }"

    input:
    tuple val(meta), path (tbl)
    
    output:
    tuple val(meta), path ("${meta.tool}.query_2_matches.json"), emit: combined
    //path ("versions.yml")   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    touch ${meta.tool}.combined.query_2_matches.json
    agg=${meta.tool}.combined.query_2_matches.json    

    for file in *${meta.tool}.query_2_matches.json; do
        combine_outputs.py \$agg \$file
    done
    mv \$agg ${meta.tool}.query_2_matches.json

    """

    stub:
    """
    touch test.txt
    """
}

