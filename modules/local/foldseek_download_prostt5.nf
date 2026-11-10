// One-time download of the ProstT5 weights foldseek uses to predict 3Di
// structural-alphabet sequences directly from amino acid sequence (no
// folding required). This is only invoked when
// ${params.prostt5_weights_dir}/.download_complete doesn't already exist
// (see the gate in workflows/orthologymap.nf) -- the weights are a multi-GB
// one-time download, not something to repeat on every run.

process FOLDSEEK_DOWNLOAD_PROSTT5 {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek' :
        'quay.io/rgrindle/foldseek:latest-cuda12' }"

    output:
    path "weights", emit: weights
    path ".download_complete", emit: marker
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tmp
    foldseek_avx2 databases ProstT5 weights tmp --threads \$(nproc)
    touch .download_complete

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek_avx2 version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p weights
    touch .download_complete
    touch versions.yml
    """
}
