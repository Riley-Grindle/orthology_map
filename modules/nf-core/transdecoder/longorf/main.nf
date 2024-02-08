process TRANSDECODER_LONGORF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
    '/comp-bio-aging/transdecoder' }"

    input:
    tuple val(meta), path(fasta)
    path(genetx_map)    

    output:
    tuple val(meta), path("*/*.pep") , emit: pep
    path "versions.yml"                       , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    TransDecoder.LongOrfs \\
        $args \\
        -t \\
        $fasta \\
        --gene_trans_map $genetx_map

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.LongOrfs --version) | sed -e "s/TransDecoder.LongOrfs //g")
    END_VERSIONS
    """
}
