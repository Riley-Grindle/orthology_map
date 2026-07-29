
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
    # Build i-ADHoRe settings file
    # ------------------------------------------------------------------
    {
        echo "genome= ${query_id}"
        for f in \$(ls ${query_genelists}/*.lst | sort); do
            echo "  \$f"
        done
        echo ""
        echo "genome= ${ref_id}"
        for f in \$(ls ${ref_genelists}/*.lst | sort); do
            echo "  \$f"
        done
        echo ""
        echo "blast_table= ${blast}"
        echo "output_path= output/"
        echo ""
        echo "number_of_anchors= 3"
        echo "gap_size= 30"
        echo "cluster_gap= 35"
        echo "q_value= 0.75"
        echo "prob_cutoff= 0.001"
        echo "anchor_points= 3"
        echo "alignment_method= gg2"
        echo "level_2_only= false"
        echo "table_type= family"
        echo "multiplicon_inheritance= true"
        echo "number_of_threads= ${task.cpus}"
    } > iadhore.ini

    # ------------------------------------------------------------------
    # Run i-ADHoRe
    # ------------------------------------------------------------------
    i-adhore iadhore.ini ${args}

    # ------------------------------------------------------------------
    # Convert anchorpoints to voting-compatible TSV.
    # anchorpoints.txt columns: id  multiplicon  gene_x  coord_x  gene_y  coord_y  is_tandem
    # gene_x = query (first genome in settings), gene_y = reference (second genome)
    # ------------------------------------------------------------------
    printf 'QUERY_ID\tREF_ID\tSCORE\n' > ${out_tsv}
    if [ -f output/anchorpoints.txt ]; then
        tail -n +2 output/anchorpoints.txt | \\
            awk 'BEGIN{OFS="\t"} NF>=5 {print \$3, \$5, 1.0}' \\
            >> ${out_tsv}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        i-adhore: \$(i-adhore --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${query_meta.id}.${ref_meta.id}.synteny_std.tsv
    touch versions.yml
    """
}
