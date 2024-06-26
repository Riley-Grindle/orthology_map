
process PRE_PROC {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/pre_proc' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path ("fasta_dir")

    script:
    """
    
    FILE_NAME="\$(basename $gtf | cut -d. -f1)"
    fasta_header_format.py $fasta $gtf "\$FILE_NAME".formatted.fasta
    
    mkdir fasta_dir
    cp *.formatted.fasta fasta_dir/
    """

}
