


process DIAMOND {
    tag "Making Diamond db using: $fasta"
    label "process_high"

    conda "conda-forge::python=3.9.5"
    // This database is only ever read by EGGNOGMAPPER's bundled diamond
    // binary (via --dmnd_db), so it must be built with a diamond binary
    // compatible with that one -- diamond's .dmnd format changed in a
    // breaking way between the 0.9.x and 2.x series, and rgrindle/diamond
    // (0.9.14) is incompatible with eggnogg-mapper's bundled 2.x diamond
    // ("Database was built with an older version of Diamond and is
    // incompatible"). Building it with the exact same image that will read
    // it guarantees compatibility and stays in sync automatically if that
    // image is ever updated.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/eggnogg-mapper' }"

    input:
    tuple val(meta), path (fasta)


    output:
    path ("*.dmnd"), emit: egg
    //path ("versions.yml")            , emit: versions

    script:
    def args = task.ext.args  ?: ''
    """
    diamond makedb --in $fasta -d ${meta.id}

    """

    stub:
    """
    touch test.txt
    """
}

