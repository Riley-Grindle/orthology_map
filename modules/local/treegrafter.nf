

process TREEGRAFTER {

    tag "Grafting $project_id using: $fasta"
    label "process_high"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/treegrafter' }"

    input:
    path fasta
    path sup_info
    val project_id

    output:
    path ("*.out.txt"), emit: tree
    //path ("versions.yml")    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    perl /opt/TreeGrafter.git/treeGrafter.pl \\
        -f $fasta \\
        -o $project_id\.out.txt \\
        -cpus ${task.cpus} \\
        -d $sup_info -auto \\
        ${args}

    """
    

    stub:
    """
    touch $fasta
    """

}
