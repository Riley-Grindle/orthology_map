// Build one foldseek structure DB per genome (query, or one reference),
// for the structural-comparison subworkflow.
//
// Two coverage sources are combined, per protein:
//   - a real structure, if the user supplied one in pdb_dir (alphafold mode
//     only) and it's named <transcript_id>.pdb
//   - otherwise, a 3Di sequence predicted directly from the amino acid
//     sequence by foldseek's built-in ProstT5 model (no folding at all) --
//     this is also the ONLY path used when fold_method == 'prostt5' (the
//     default), in which case pdb_dir is never consulted.
//
// foldseek's own createdb handles both structure input (a directory of PDB
// files) and sequence input (a FASTA + --prostt5-model) with the same
// command; we just need to merge the two resulting DBs into one when a
// genome has both kinds of coverage. A structure DB is actually several
// parallel sequence-DB "objects" sharing one id space: the AA sequences
// themselves, a header DB (_h), the predicted/derived 3Di sequences (_ss),
// and (PDB-derived entries only) real C-alpha coordinates (_ca). We merge
// the first three with foldseek's own `concatdbs` (which correctly
// renumbers/reindexes them) and deliberately DROP _ca when mixing: without
// it foldseek just can't compute coordinate-dependent extras (TM-score,
// "prob", LDDT) for the combined DB, but the core 3Di/AA alignment search
// doesn't need it at all. FOLDSEEK_SEARCH's output columns are chosen
// accordingly (see that module).
//
// concatdbs itself doesn't touch the plain-text .lookup/.source files
// (they're not index-based objects), so those are merged by hand -- the
// second DB's internal ids need to be offset by the first DB's entry count
// before concatenating, since createdb always numbers each DB from 0.

process FOLDSEEK_CREATEDB {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/foldseek' :
        'quay.io/rgrindle/foldseek:latest-cuda12' }"

    input:
    tuple val(meta), path(pep_fasta)
    path  prostt5_weights
    path  pdb_dir
    val   fold_method

    output:
    tuple val(meta), path("${meta.id}.structDB*"), emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export PATH="${projectDir}/bin:\$PATH"

    merge_object() {
        # \$1 = object suffix, e.g. "", "_h", "_ss"
        local suffix="\$1"
        if [ -f "pdb_db\${suffix}" ] && [ -f "seq_db\${suffix}" ]; then
            foldseek_avx2 concatdbs "pdb_db\${suffix}" "seq_db\${suffix}" "${meta.id}.structDB\${suffix}" --threads \$(nproc)
        elif [ -f "pdb_db\${suffix}" ]; then
            for f in "pdb_db\${suffix}" "pdb_db\${suffix}.index" "pdb_db\${suffix}.dbtype"; do
                [ -f "\$f" ] && cp "\$f" "${meta.id}.structDB\${suffix}\${f#pdb_db\${suffix}}"
            done
        else
            for f in "seq_db\${suffix}" "seq_db\${suffix}.index" "seq_db\${suffix}.dbtype"; do
                [ -f "\$f" ] && cp "\$f" "${meta.id}.structDB\${suffix}\${f#seq_db\${suffix}}"
            done
        fi
    }

    HAVE_PDB=0
    HAVE_SEQ=0

    if [ "${fold_method}" = "alphafold" ] && [ -d "${pdb_dir}" ]; then
        mkdir -p has_pdb
        split_pdb_coverage.py ${pep_fasta} ${pdb_dir} has_pdb missing.fasta
        if [ -n "\$(ls -A has_pdb 2>/dev/null)" ]; then
            foldseek_avx2 createdb has_pdb pdb_db --threads \$(nproc) ${args}
            HAVE_PDB=1
        fi
        if [ -s missing.fasta ]; then
            foldseek_avx2 createdb missing.fasta seq_db --prostt5-model ${prostt5_weights} --gpu 1 --threads \$(nproc) ${args}
            HAVE_SEQ=1
        fi
    else
        foldseek_avx2 createdb ${pep_fasta} seq_db --prostt5-model ${prostt5_weights} --gpu 1 --threads \$(nproc) ${args}
        HAVE_SEQ=1
    fi

    merge_object ""
    merge_object "_h"
    merge_object "_ss"

    if [ "\$HAVE_PDB" = "1" ] && [ "\$HAVE_SEQ" = "1" ]; then
        offset=\$(wc -l < pdb_db.lookup)
        cat pdb_db.lookup > ${meta.id}.structDB.lookup
        awk -v off="\$offset" 'BEGIN{OFS="\t"} {\$1=\$1+off; print}' seq_db.lookup >> ${meta.id}.structDB.lookup
        cat pdb_db.source seq_db.source > ${meta.id}.structDB.source
    elif [ "\$HAVE_PDB" = "1" ]; then
        cp pdb_db.lookup ${meta.id}.structDB.lookup
        cp pdb_db.source ${meta.id}.structDB.source
    else
        cp seq_db.lookup ${meta.id}.structDB.lookup
        cp seq_db.source ${meta.id}.structDB.source
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        foldseek: \$(foldseek_avx2 version)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.structDB
    touch versions.yml
    """
}
