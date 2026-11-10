// Structural comparison: query genome's foldseek DB vs one reference
// genome's foldseek DB. Standalone diagnostic output, not fed into
// CONSENSUS_VOTE (structural-similarity signal quality hasn't been
// validated against the other evidence sources yet).
//
// Output columns are deliberately limited to ones foldseek can compute from
// the 3Di/AA alignment alone: when a genome's DB mixes PDB-derived and
// ProstT5-derived entries (alphafold mode with partial coverage) there are
// no real coordinates for the ProstT5 entries, and columns that need them
// (prob, lddt, alntmscore, rmsd, ...) make convertalis fail outright rather
// than just omitting them for those rows.

process FOLDSEEK_SEARCH {
    tag "${query_meta.id}:${ref_meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek' :
        'quay.io/rgrindle/foldseek:latest-cuda12' }"

    input:
    tuple val(query_meta), path(query_db)
    tuple val(ref_meta),   path(ref_db)

    output:
    tuple val(ref_meta), path("*.structural_hits.tsv"), emit: hits
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def out_tsv = "${query_meta.id}.${ref_meta.id}.structural_hits.tsv"
    """
    mkdir -p tmp
    foldseek_avx2 search ${query_meta.id}.structDB ${ref_meta.id}.structDB result tmp --threads \$(nproc) ${args}
    foldseek_avx2 convertalis ${query_meta.id}.structDB ${ref_meta.id}.structDB result ${out_tsv} \\
        --format-mode 4 \\
        --format-output query,target,evalue,bits,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend \\
        --threads \$(nproc)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek_avx2 version)
    END_VERSIONS
    """

    stub:
    """
    touch ${query_meta.id}.${ref_meta.id}.structural_hits.tsv
    touch versions.yml
    """
}
