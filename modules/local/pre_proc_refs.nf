params.refs = ""


process PRE_PROC {
    tag "Formatting Uniprot Reference Headers"
    label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/pre_proc' }"

    input:
    path references

    output:
    path ("fasta_dir")

    script:
    """
    for file in "$references"/*; do
        FILE_NAME="\$(basename "\$file" | cut -d. -f1)"
        fasta_header_format.py \$file "\$FILE_NAME".formatted.fasta
    done
    mkdir fasta_dir
    cp *.formatted.fasta fasta_dir/
    """

}


workflow {

    ch_refs = PRE_PROC{ params.refs } 
    ch_refs.view()

}
