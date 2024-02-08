
process GTF_2_GENETX_MAP {
    tag "Mapping input transcripts to gene names"
    label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/gtf_2_genetx_map' }"

    input:
    path gtf

    output:
    path ("output_mapping.txt"), emit: genetx_map

    script:
    """
    gtf2gene-tx-map.R $gtf
    """

}
