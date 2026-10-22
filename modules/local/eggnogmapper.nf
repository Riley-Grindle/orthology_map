


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
    path dmnd

    output:
    path ("*.emapper.seed_orthologs"), emit: egg
    path ("versions.yml")            , emit: versions

    script:
    def args = task.ext.args  ?: ''
    if (dmnd instanceof java.util.ArrayList){
        """
        emapper.py \\
            -i $fasta/ \\
            -m diamond \\
            -o $project_id \\
            --itype proteins \\
            --cpu \$(nproc) \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            eggnog_mapper: \$(emapper.py --version)
        END_VERSIONS
        """
    } else {
        // `-d`/`--database` selects eggnog-mapper's HMMER-mode database, not
        // a custom diamond database -- with `-m diamond` it's silently
        // ignored and emapper.py falls back to its bundled, full
        // cross-species eggnog_proteins.dmnd. `--dmnd_db` is the diamond-mode
        // flag for pointing at this reference-species-only database.
        """
        emapper.py \\
            -i $fasta/ \\
            -m diamond \\
            -o $project_id \\
            --itype proteins \\
            --cpu \$(nproc) \\
            --dmnd_db $dmnd \\
            ${args}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            eggnog_mapper: \$(emapper.py --version)
        END_VERSIONS
        """
    }

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog_mapper: \$(emapper.py --version)
    END_VERSIONS
    """
}

