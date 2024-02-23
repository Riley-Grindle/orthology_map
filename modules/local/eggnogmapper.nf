


process EGGNOGMAPPER {
    tag "Eggnogging $project_id using: $fasta"
    label "process_high"    
    
    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/eggnogg-mapper' }"

    input:
    path fasta
    val project_id

    output:
    path ("*.emapper.seed_orthologs"), emit: egg
    path ("versions.yml")            , emit: versions

    script:
    def args = task.ext.args  ?: ''
    """
    emapper.py \\
        -i $fasta/ \\
        -m diamond \\
        -o $project_id \\
        --itype proteins \\
        --cpu ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog_mapper: \$(emapper.py --version)
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog_mapper: \$(emapper.py --version)
    END_VERSIONS
    """
}

