


process PREP_INPUT {
    tag "Staging ${meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(ref_fastas)
    tuple val(meta_1), path(query_fasta)
    val project_id

    output:
    tuple val(meta), path ("./ortho_f"), emit: ortho_f
    tuple val(meta), path ("./ortho_l/odbwork"), emit: ortho_l_work
    tuple val(meta), path ("./ortho_l/odbdata"), emit: ortho_l_data
    path "./eggnog/*", emit: egg
    path "./treegrafter/*", emit: tree
    path "blank.txt", emit: blank

    
    script:
    def query_fasta_ext = "$query_fasta" + ".fasta"
    """
    mkdir ortho_f; mkdir ortho_l; mkdir treegrafter; mkdir eggnog
    mv $query_fasta $query_fasta_ext
    cd ortho_l 
    mkdir odbwork
    mkdir odbdata
    cd ..
    touch blank.txt
    mv blank.txt ./ortho_l/odbwork/
    cp -r $ref_fastas/* ./ortho_f/; cp $query_fasta_ext ./ortho_f/
    cp -r $ref_fastas/* ./ortho_l/odbdata/; cp $query_fasta_ext ./ortho_l/odbdata/
    cp $query_fasta_ext ./treegrafter/ 
    cp $query_fasta_ext ./eggnog/
    touch blank.txt
    """

    stub:
    """
    touch $fastas
    """
}

