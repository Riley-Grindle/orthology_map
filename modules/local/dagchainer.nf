

process DAGCHAINER {
    
    tag "Synteny Identification: $project_id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/dagchainer:latest' }"

    input:
    path(blast_tbl)
    path(query_gtf)
    path(ref_gtf)
    val(project_id)

    output:
    path ("*.aligncoords")               , emit: synteny
    path ("versions.yml")                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    blast_2_dag.py $query_gtf $ref_gtf $blast_tbl    

    filter_repetitive_matches.pl ${args} < dagchainer.db.tsv > dagchainer.db.filtered.tsv
    
    run_DAG_chainer.pl -i dagchainer.db.filtered.tsv -s -I


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
