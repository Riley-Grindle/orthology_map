
// All-vs-all protein similarity search between the query and one reference species.
// Output is a 3-column blast table (qseqid, sseqid, evalue) compatible with i-ADHoRe.

process DIAMOND_BLASTP {
    tag "${query_meta.id}:${ref_meta.id}"
    label 'process_high'

    conda "bioconda::diamond=2.1.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.8--h43eeafb_0' :
        'quay.io/biocontainers/diamond:2.1.8--h43eeafb_0' }"

    input:
    tuple val(query_meta), path(query_pep)
    tuple val(ref_meta),   path(ref_pep)

    output:
    tuple val(ref_meta), path("${query_meta.id}.${ref_meta.id}.blast"), emit: blast
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def out_file = "${query_meta.id}.${ref_meta.id}.blast"
    """
    diamond makedb --in ${ref_pep} -d ref_db --quiet

    diamond blastp \\
        -q ${query_pep} \\
        -d ref_db \\
        -o ${out_file} \\
        --outfmt 6 qseqid sseqid evalue \\
        -p ${task.cpus} \\
        -e 1e-5 \\
        --more-sensitive \\
        --quiet \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | head -1 | sed 's/diamond version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${query_meta.id}.${ref_meta.id}.blast
    touch versions.yml
    """
}
