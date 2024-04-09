
process TRANSDECODER_LONGORF {
    
    tag '${meta.id}'
    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4' :
    'rgrindle/transdecoder' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta_1), path(genetx_map)    

    output:
    tuple val(meta), path("*transdecoder_dir/*.pep")  , emit: pep
    tuple val(meta), path("*.transcripts.fa")         , emit: transcripts
    tuple val(meta), path("*transdecoder_dir")        , emit: out_dir
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tx_fasta = "${meta.id}.transcripts.fa"
    """
    cp $fasta ${meta.id}.transcripts.fa

    TransDecoder.LongOrfs \\
        $args \\
        -t \\
        $tx_fasta \\
        --gene_trans_map $genetx_map
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.LongOrfs --version) | sed -e "s/TransDecoder.LongOrfs //g")
    END_VERSIONS
    """
}
