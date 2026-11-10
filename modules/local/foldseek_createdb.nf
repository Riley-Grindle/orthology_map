// Build one foldseek structure DB per genome (query, or one reference),
// for the structural-comparison subworkflow.
//
// Up to three coverage sources are combined, per protein, in priority
// order:
//   1. a real structure, if the user supplied one in pdb_dir (alphafold
//      mode only) and it's named <transcript_id>.pdb
//   2. a previously-computed ProstT5 3Di prediction, if one already exists
//      in the persistent prostt5_3di_cache_dir (any mode) -- reused
//      without ever invoking ProstT5 again
//   3. otherwise, a fresh 3Di sequence predicted directly from the amino
//      acid sequence by foldseek's built-in ProstT5 model -- the only
//      path used when fold_method == 'prostt5' (the default) and nothing
//      is cached yet. Newly-predicted 3Di sequences are unpacked and
//      published back into prostt5_3di_cache_dir (see conf/modules.config)
//      so future runs -- any species, any project -- never re-predict
//      that exact protein again.
//
// foldseek's own createdb handles both structure input (a directory of PDB
// files) and sequence input (a FASTA + --prostt5-model) with the same
// command; we just need to merge however many of the three DBs apply into
// one. A structure DB is actually several parallel sequence-DB "objects"
// sharing one id space: the AA sequences themselves, a header DB (_h), the
// predicted/derived 3Di sequences (_ss), and (PDB-derived entries only)
// real C-alpha coordinates (_ca). We merge the first three with foldseek's
// own `concatdbs` (which correctly renumbers/reindexes them) and
// deliberately DROP _ca when mixing: without it foldseek just can't
// compute coordinate-dependent extras (TM-score, "prob", LDDT) for the
// combined DB, but the core 3Di/AA alignment search doesn't need it at
// all. FOLDSEEK_SEARCH's output columns are chosen accordingly (see that
// module).
//
// The cache tier works because a sequence-DB object's on-disk format
// doesn't care what alphabet its text is in -- confirmed empirically, a
// real "_ss" object has the exact same dbtype bytes (0 0 0 0) and entry
// format (raw text + "\n\0" per entry) as foldseek's ordinary AA sequence
// objects. BUT foldseek's own `createdb`, called without --prostt5-model
// and without real structure input, does NOT build a plain sequence DB at
// all -- it always takes the structure-parsing path and exits non-zero
// ("No structures found in given input") with empty output, even on a
// perfectly good FASTA. So bin/fasta_to_seqdb.py builds the handful of
// flat files (data/.index/.dbtype/.lookup/.source, plus the _h header
// companion) directly, bypassing createdb entirely for the cache tier.
// The matching AA/header objects for those same proteins are built the
// same way from their original amino acid sequences, in lockstep order
// with the 3Di fasta (see bin/split_3di_cache.py) so the two DBs' internal
// ids line up 1:1 before merging.
//
// concatdbs itself doesn't touch the plain-text .lookup/.source files
// (they're not index-based objects), so those are merged by hand -- each
// later source's internal ids need to be offset by the running entry count
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
    path  prostt5_3di_cache_dir
    path  pdb_dir
    val   fold_method

    output:
    tuple val(meta), path("${meta.id}.structDB*"), emit: db
    path "*.3di", emit: new_cache_entries, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export PATH="${projectDir}/bin:\$PATH"

    merge_lookup_source() {
        # Append, offsetting internal ids by the running entry count, for
        # whichever of pdb_db/cached_db/fresh_db exist.
        : > "${meta.id}.structDB.lookup"
        : > "${meta.id}.structDB.source"
        offset=0
        for src in pdb_db cached_db fresh_db; do
            if [ -f "\${src}.lookup" ]; then
                awk -v off="\$offset" 'BEGIN{OFS="\t"} {\$1=\$1+off; print}' "\${src}.lookup" >> "${meta.id}.structDB.lookup"
                cat "\${src}.source" >> "${meta.id}.structDB.source"
                offset=\$((offset + \$(wc -l < "\${src}.lookup")))
            fi
        done
    }

    merge_object() {
        # \$1 = object suffix, e.g. "", "_h", "_ss"
        suffix="\$1"
        srcs=""
        for src in pdb_db cached_db fresh_db; do
            [ -f "\${src}\${suffix}" ] && srcs="\$srcs \$src"
        done
        set -- \$srcs
        if [ \$# -eq 0 ]; then
            return
        elif [ \$# -eq 1 ]; then
            for f in "\$1\${suffix}" "\$1\${suffix}.index" "\$1\${suffix}.dbtype"; do
                [ -f "\$f" ] && cp "\$f" "${meta.id}.structDB\${suffix}\${f#\$1\${suffix}}"
            done
        else
            acc="\$1\${suffix}"
            shift
            for src in "\$@"; do
                tmp="merged_\${src}\${suffix}"
                foldseek_avx2 concatdbs "\$acc" "\${src}\${suffix}" "\$tmp" --threads \$(nproc)
                acc="\$tmp"
            done
            for f in "\$acc" "\$acc.index" "\$acc.dbtype"; do
                [ -f "\$f" ] && cp "\$f" "${meta.id}.structDB\${suffix}\${f#\$acc}"
            done
        fi
    }

    # ------------------------------------------------------------------
    # Tier 1: user-supplied real PDB structures (alphafold mode only)
    # ------------------------------------------------------------------
    REMAINING_FASTA="${pep_fasta}"
    if [ "${fold_method}" = "alphafold" ] && [ -d "${pdb_dir}" ]; then
        mkdir -p has_pdb
        split_pdb_coverage.py ${pep_fasta} ${pdb_dir} has_pdb missing_after_pdb.fasta
        if [ -n "\$(ls -A has_pdb 2>/dev/null)" ]; then
            foldseek_avx2 createdb has_pdb pdb_db --threads \$(nproc) ${args}
        fi
        REMAINING_FASTA=missing_after_pdb.fasta
    fi

    # ------------------------------------------------------------------
    # Tier 2: previously-computed ProstT5 3Di predictions (any mode)
    # ------------------------------------------------------------------
    if [ -d "${prostt5_3di_cache_dir}" ] && [ -s "\$REMAINING_FASTA" ]; then
        split_3di_cache.py "\$REMAINING_FASTA" ${prostt5_3di_cache_dir} cached_aa.fasta cached_3di.fasta needs_prostt5.fasta
        if [ -s cached_aa.fasta ]; then
            fasta_to_seqdb.py cached_aa.fasta cached_db
            fasta_to_seqdb.py cached_3di.fasta cached_db_ss_raw
            for f in cached_db_ss_raw cached_db_ss_raw.index cached_db_ss_raw.dbtype; do
                [ -f "\$f" ] && mv "\$f" "cached_db_ss\${f#cached_db_ss_raw}"
            done
        fi
    else
        cp "\$REMAINING_FASTA" needs_prostt5.fasta
    fi

    # ------------------------------------------------------------------
    # Tier 3: fresh ProstT5 prediction for anything still uncovered
    # ------------------------------------------------------------------
    if [ -s needs_prostt5.fasta ]; then
        foldseek_avx2 createdb needs_prostt5.fasta fresh_db --prostt5-model ${prostt5_weights} --gpu 1 --threads \$(nproc) ${args}

        # Publish the newly-predicted 3Di sequences back into the
        # persistent cache so future runs skip ProstT5 for these proteins.
        # Unpacked flat into the task's own root (not a subdirectory) so
        # publishDir's pattern = '*.3di' lands them flat in
        # prostt5_3di_cache_dir too, matching the plain <accession>.3di
        # layout split_3di_cache.py looks up.
        cp fresh_db.lookup fresh_db_ss.lookup
        foldseek_avx2 unpackdb fresh_db_ss . --unpack-name-mode 1 --unpack-suffix .3di --threads \$(nproc)
    fi

    merge_object ""
    merge_object "_h"
    merge_object "_ss"
    merge_lookup_source

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
