

process VOTE_BEST_MATCH {
    
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/tool_vote:latest' }"

    input:
    path jsons
    path query_fasta

    output:
    path ("voted_orthologs.json"), emit: combined_orthologs
    //path ("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    touch voted_orthologs.json
    grep ">" $query_fasta > queries.txt

    for file in *query_2_matches.json; do
        vote.py queries.txt voted_orthologs.json \$file
    done
   
    species_layer.py voted_orthologs.json    
    """

    stub:
    """
    
    """
}

