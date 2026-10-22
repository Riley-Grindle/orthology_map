
// Standardize i-ADHoRe's raw synteny_std.tsv (QUERY_ID, REF_ID, SCORE) into
// the <QUERY, MATCH, GENE, SPECIES> format CONSENSUS_VOTE:TBL_2_JSON expects,
// the same way POST_PROC does for the ortho_f/ortho_l/eggnog/tree evidence
// sources. Without this, TBL_2_JSON's "expected >=4 tab-separated fields"
// check fails on synteny's raw 3-column output.

process SYNTENY_POST {
    tag "${ref_meta.id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'rgrindle/post_proc' }"

    input:
    tuple val(ref_meta), path(synteny_tsv)
    path ensembl_data    // absolute path or null — not staged, accessed directly on shared fs

    output:
    tuple val(ref_meta), path("*.synteny_std.tsv"), emit: synteny

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    if [ -n "$ensembl_data" ]; then
        synteny_post.py ${synteny_tsv} $ensembl_data /sup_data/prefix_2_file.json /sup_data/prefix_2_species.json ${ref_meta.id}.synteny_std.tsv
    else
        awk 'BEGIN{FS=OFS="\t"} NR==1{print "QUERY","MATCH","GENE","SPECIES"; next} {print \$1,\$2,"NA","NA"}' ${synteny_tsv} > ${ref_meta.id}.synteny_std.tsv
    fi
    """
}
