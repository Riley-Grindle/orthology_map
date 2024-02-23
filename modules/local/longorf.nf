
process TRANSDECODER_LONGORF {

    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
    'rgrindle/transdecoder' }"

    input:
    path(fasta)
    path(genetx_map)    

    output:
    path("*transdecoder_dir/*.pep")  , emit: pep
    path("*transdecoder_dir")        , emit: out_dir
    path "versions.yml"              , emit: versions

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
