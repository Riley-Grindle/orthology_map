
process PANTHER_API {
    tag "Pulling gene ids from PantherDB"
    label "process_medium"

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/panther_api_sf' }"

    input:
    path tree_out

    output:
    path "*orthologs.txt"

    script:
    """
    bash /post_panther/post_panther.sh
    mv orthologs.txt "\$(basename $tree_out)_orthologs.txt"
    """

}
