


process ORTHOLOGER {

    tag "Logging $project_id"
    label 'process_medium'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/orthologer' }"
    

    input:
    tuple val(meta), path (odbdata)
    tuple val(meta), path (odbwork)
    val path_to_work
    val project_id 

    output:
    tuple val(meta), path ("$odbwork/Cluster/odbproj.og"), emit: loger
    path ("$odbwork/versions.yml")               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args  ?: ''
    """
    WORKDIR=\$(pwd)

    # ORTHOLOGER derives its internal species labels from the part of each
    # input filename before the first '.', and only allows [a-z0-9_]
    # (case-insensitive). Sanitize filenames here so labels like
    # "\${params.project_id}" (which may contain '-' or '.') don't trip its
    # "bad characters" check, without changing the actual pipeline parameter.
    mkdir -p /odbdata
    for f in $odbdata/*; do
        fname=\$(basename "\$f")
        base="\${fname%%.*}"
        ext="\${fname#*.}"
        safe_base=\$(printf '%s' "\$base" | tr -c 'A-Za-z0-9_' '_')
        if [ "\$ext" = "\$fname" ]; then
            cp -r "\$f" "/odbdata/\${safe_base}"
        else
            cp -r "\$f" "/odbdata/\${safe_base}.\${ext}"
        fi
    done
    cd /odbwork
    setup_odb.sh
    orthologer -c import
    cd Sequences/
    capitalize_files.py
    cd ../
    orthologer \\
        -c run \\
        ${args}

    # Copy results back into the Nextflow work directory so outputs are properly staged
    mkdir -p \${WORKDIR}/${odbwork}/Cluster
    cp -r Cluster/. \${WORKDIR}/${odbwork}/Cluster/

    # `orthologer.sh -v` prints a multi-line startup banner rather than a
    # clean version string; pull just the "pipeline version" line so the
    # versions.yml we write stays valid single-line YAML.
    ORTHOLOGER_VERSION=\$(./orthologer.sh -v 2>&1 | grep -m1 'pipeline version' | sed -E 's/^pipeline version[[:space:]]*:[[:space:]]*//' || true)
    ORTHOLOGER_VERSION=\${ORTHOLOGER_VERSION:-unknown}
cat <<-END_VERSIONS > \${WORKDIR}/${odbwork}/versions.yml
"${task.process}":
    orthologer: "\${ORTHOLOGER_VERSION}"
END_VERSIONS
    """

    stub:
    """
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    orthologer: "unknown"
END_VERSIONS
    """

}
