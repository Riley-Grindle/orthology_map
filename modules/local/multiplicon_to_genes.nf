// Expand IADHORE's multiplicons.txt into a gene-level diagnostic table.
// Split out from IADHORE itself because rgrindle/i-adhore:latest (Ubuntu
// 20.04, i-ADHoRe's own build image) has no python interpreter at all --
// only perl/awk -- so this needs a container that actually has python3.

process MULTIPLICON_TO_GENES {
    tag "${ref_meta.id}"
    label 'process_low'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/eggnogg-mapper' }"

    input:
    tuple val(ref_meta), path(iadhore_out), path(ref_genelists)
    tuple val(query_meta), path(query_genelists)

    output:
    path "multiplicon_genes.txt", emit: genes, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    if [ -f ${iadhore_out}/multiplicons.txt ]; then
        multiplicon_to_genes.py ${iadhore_out} multiplicon_genes.txt \\
            ${query_meta.id} ${query_genelists} \\
            ${ref_meta.id} ${ref_genelists}
    fi
    """

    stub:
    """
    touch multiplicon_genes.txt
    """
}
