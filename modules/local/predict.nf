
process TRANSDECODER_PREDICT {
    tag '${meta.id}'
    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4':
        'rgrindle/transdecoder' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(pfam)
    tuple val(meta), path(long_orf_dir)
    val(project_id)

    output:
    path("*.initial.transdecoder.pep") , emit: orig
    tuple val(meta), path("*.final.transdecoder.pep")  , emit: pep
    path("*.transdecoder.gff3") , emit: gff3
    path("*.transdecoder.cds")  , emit: cds
    path("*.transdecoder.bed")  , emit: bed
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    TransDecoder.Predict \\
        $args \\
        -t \\
        $fasta \\
        --retain_pfam_hits $pfam \\
        ${args}

    mv *.transdecoder.pep ${meta.id}.initial.transdecoder.pep
    sed 's/>[^ ]*//' *.initial.transdecoder.pep | sed 's/^ />/'  > ${meta.id}.final.transdecoder.pep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.Predict --version) | sed -e "s/TransDecoder.Predict //g")
    END_VERSIONS
    
    """
}
