
// i-ADHoRe pairwise synteny between the query species and one reference species.
// Runs once per reference in the samplesheet via channel broadcasting.
//
// Inputs:
//   query_genelists — directory of per-scaffold .lst files (transcript_id+/-)
//   ref_genelists   — same, for the reference species
//   blast           — 3-column DIAMOND blast table (qseqid sseqid evalue)
//
// Output synteny_std.tsv columns: QUERY_ID  REF_ID  SCORE
//   where QUERY_ID/REF_ID are transcript_ids matching TransDecoder final.pep headers.

process IADHORE {
    tag "${query_meta.id}:${ref_meta.id}"
    label 'process_high'

    conda "bioconda::i-adhore=3.0.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'rgrindle/i-adhore:latest' :
        'rgrindle/i-adhore:latest' }"

    input:
    tuple val(query_meta), path(query_genelists)
    tuple val(ref_meta),   path(ref_genelists), path(blast)

    output:
    tuple val(ref_meta), path("*.synteny_std.tsv"), emit: synteny
    tuple val(ref_meta), path("output")            , emit: outdir
    path "output/anchorpoints.txt"                 , emit: anchorpoints, optional: true
    path "output/multiplicons.txt"                 , emit: multiplicons, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def query_id = query_meta.id
    def ref_id   = ref_meta.id
    def out_tsv  = "${query_id}.${ref_id}.synteny_std.tsv"
    """
    mkdir -p output

    # ------------------------------------------------------------------
    # i-ADHoRe's "pairs" table reader only splits each line on the FIRST
    # tab (geneA = up to tab 1, geneB = rest of line), so a 3rd (evalue)
    # column gets appended onto geneB and breaks every match -- keep only
    # the first two columns. It also matches blast IDs against the
    # gene-list (.lst) IDs by exact string equality: the gene lists use
    # "gene_id~~transcript_id" (matching TransDecoder's --gene_trans_map
    # headers) with no ORF suffix, but the blast table's IDs carry
    # TransDecoder's per-ORF ".pN" suffix (e.g. "...~~ENSDART...159919.p1").
    # Strip that too so both sides actually match.
    # ------------------------------------------------------------------
    cut -f1,2 ${blast} | sed -E 's/\\.p[0-9]+\t/\t/; s/\\.p[0-9]+\$//' > blast_table.stripped.tsv

    # ------------------------------------------------------------------
    # Build i-ADHoRe settings file
    # ------------------------------------------------------------------
    {
        echo "genome= ${query_id}"
        for f in \$(ls ${query_genelists}/*.lst | sort); do
            echo "  \$(basename "\$f" .lst) \$f"
        done
        echo ""
        echo "genome= ${ref_id}"
        for f in \$(ls ${ref_genelists}/*.lst | sort); do
            echo "  \$(basename "\$f" .lst) \$f"
        done
        echo ""
        echo "blast_table= blast_table.stripped.tsv"
        echo "output_path= output/"
        echo ""
        echo "gap_size= 30"
        echo "cluster_gap= 35"
        echo "q_value= 0.75"
        echo "prob_cutoff= 0.001"
        echo "anchor_points= 3"
        echo "alignment_method= gg2"
        echo "level_2_only= false"
        echo "table_type= pairs"
        echo "number_of_threads= \$(nproc)"
    } > iadhore.ini

    # ------------------------------------------------------------------
    # Run i-ADHoRe
    # ------------------------------------------------------------------
    i-adhore iadhore.ini ${args}

    # ------------------------------------------------------------------
    # Convert anchorpoints to voting-compatible TSV.
    # anchorpoints.txt columns: id  multiplicon  basecluster  gene_x  gene_y  coord_x  coord_y  is_real_anchorpoint
    # gene_x/gene_y are NOT consistently query/reference -- which side is
    # which flips per multiplicon -- so resolve orientation per row against
    # the query gene-list IDs rather than assuming a fixed column order.
    # ------------------------------------------------------------------
    printf 'QUERY_ID\tREF_ID\tSCORE\n' > ${out_tsv}
    if [ -f output/anchorpoints.txt ]; then
        awk 'BEGIN{OFS="\t"}
            FILENAME==ARGV[ARGC-1] {
                if (FNR==1) next
                if (NF<7) next
                gx=\$4; gy=\$5
                if (gx in q) print gx, gy, 1.0
                else if (gy in q) print gy, gx, 1.0
                next
            }
            { sub(/[+-]\$/,""); q[\$1]=1 }
        ' ${query_genelists}/*.lst output/anchorpoints.txt >> ${out_tsv}
    fi

    IADHORE_VERSION=\$(i-adhore --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1 || true)
    IADHORE_VERSION=\${IADHORE_VERSION:-unknown}
cat <<-END_VERSIONS > versions.yml
"${task.process}":
    i-adhore: "\${IADHORE_VERSION}"
END_VERSIONS
    """

    stub:
    """
    mkdir -p output
    touch ${query_meta.id}.${ref_meta.id}.synteny_std.tsv
    touch versions.yml
    """
}
