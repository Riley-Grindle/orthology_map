
process TRANSDECODER_PREDICT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::transdecoder=5.5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transdecoder:5.5.0--pl5262hdfd78af_4':
        'rgrindle/transdecoder' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(pfam)       // *.domtblout from HMMER_HMMSCAN, or [] when pfam DB not supplied
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
    def args      = task.ext.args ?: ''
    def pfam_arg  = pfam instanceof List ? '' : "--retain_pfam_hits $pfam"

    """
    # long_orf_dir is staged as a symlink into LONGORF's cached work directory.
    # TransDecoder.Predict writes its own resume checkpoints into that same
    # directory, so a stale checkpoint from a previous Predict attempt would
    # leak in here and cause this run to skip steps that never actually ran
    # in this task's work dir. Materialize a private copy and strip any
    # inherited Predict checkpoints so this run always starts clean.
    if [ -L $long_orf_dir ]; then
        real_dir="\$(readlink -f $long_orf_dir)"
        rm -f $long_orf_dir
        cp -r "\$real_dir" $long_orf_dir
    fi
    rm -rf $long_orf_dir/__checkpoints_TDpredict

    TransDecoder.Predict \\
        $args \\
        -t \\
        $fasta \\
        $pfam_arg

    mv *.transdecoder.pep ${meta.id}.initial.transdecoder.pep
    sed 's/>[^ ]*//' *.initial.transdecoder.pep | sed 's/^ />/'  > ${meta.id}.final.transdecoder.pep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(echo \$(TransDecoder.Predict --version) | sed -e "s/TransDecoder.Predict //g")
    END_VERSIONS
    """
}
