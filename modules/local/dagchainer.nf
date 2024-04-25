

process DAGCHAINER {
    
    tag "Synteny Identification: $project_id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/dagchainer:latest' }"

    input:
    path(blast_tbl)
    tuple val(meta), path(query_gtf)
    tuple val(meta), path(ref_gtf)
    val(project_id)

    output:
    path ("*.aligncoords")               , emit: synteny
    path ("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    blast_2_dag.py $query_gtf $ref_gtf $blast_tbl    

    filter_repetitive_matches.pl ${args} < dagchainer.db.tsv > dagchainer.db.filtered.tsv
    
    run_DAG_chainer.pl -i dagchainer.db.filtered.tsv -s -I ${args2}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    
    END_VERSIONS
    """
}

