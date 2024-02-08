process HMMER_HMMSCAN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::hmmer=3.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:v3.2.1dfsg-1-deb_cv1' :
        'biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1' }"

    input:
    tuple val(meta), path(hmmfile)
    path(longest_orfs)

    output:
    tuple val(meta), path(*.domtblout)
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "${prefix}.domtblout"
    """
    hmmscan \\
        $args \\
        --cpu $task.cpus \\
        --domtblout $output \\
        $hmmfile \\
        $longest_orfs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmscan -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
