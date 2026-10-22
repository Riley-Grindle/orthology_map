
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def summary_params = paramsSummaryMap(workflow)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { ORTHOFINDER } from '../modules/local/orthofinder.nf'
include { EGGNOGMAPPER } from '../modules/local/eggnogmapper.nf'
include { DIAMOND } from '../modules/local/diamond.nf'
include { TREEGRAFTER } from "../modules/local/treegrafter.nf"
include { ORTHOLOGER } from "../modules/local/orthologer.nf"
include { PRE_PROC } from "../modules/local/pre_proc_refs.nf"
include { PREP_INPUT } from "../modules/local/prep_input.nf"
include { POST_PROC } from "../modules/local/post_proc.nf"
include { STAGE_OUTS } from "../modules/local/stage_outs.nf"
include { PANTHER_API } from "../modules/local/panther_api.nf"
include { COLLECT_CHUNKS } from "../modules/local/collect_chunks.nf"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "../modules/nf-core/custom/dumpsoftwareversions/main"
include { GUNZIP as GUNZIP_REF } from "../modules/local/gunzip.nf"
include { GUNZIP as GUNZIP_IN } from "../modules/local/gunzip.nf"
include { GFFREAD } from "../modules/local/gffread/main"
include { GFFREAD as GFFREAD_1 } from "../modules/local/gffread/main"
include { GTF_TO_GENELISTS as GTF_TO_GENELISTS_QUERY } from "../modules/local/gtf_to_genelists.nf"
include { GTF_TO_GENELISTS as GTF_TO_GENELISTS_REF   } from "../modules/local/gtf_to_genelists.nf"
include { DIAMOND_BLASTP } from "../modules/local/diamond_blastp.nf"
include { IADHORE } from "../modules/local/iadhore.nf"
include { MULTIPLICON_TO_GENES } from "../modules/local/multiplicon_to_genes.nf"
include { SYNTENY_POST } from "../modules/local/synteny_post.nf"
include { FOLDSEEK_DOWNLOAD_PROSTT5 } from "../modules/local/foldseek_download_prostt5.nf"
include { FOLDSEEK_CREATEDB as FOLDSEEK_CREATEDB_QUERY } from "../modules/local/foldseek_createdb.nf"
include { FOLDSEEK_CREATEDB as FOLDSEEK_CREATEDB_REF   } from "../modules/local/foldseek_createdb.nf"
include { FOLDSEEK_SEARCH } from "../modules/local/foldseek_search.nf"

include { TRANSDECODER } from "../subworkflows/local/transdecoder.nf"
include { TRANSDECODER as REF_TRANSDECODER } from "../subworkflows/local/transdecoder.nf"
include { CONSENSUS_VOTE } from "../subworkflows/local/consensus_vote.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ORTHOLOGYMAP {

    ch_versions  = Channel.empty()

    ch_refs      = Channel
                    .fromSamplesheet("sample_sheet")

    // Map each reference species' NCBI taxon ID to the gene-ID prefix PANTHER
    // uses for that species in its family trees (e.g. 10090 -> "MGI"), so
    // PANTHER_API can restrict "genes under the graft point" to just the
    // reference species instead of every species in that part of the tree.
    // Species not in this table (e.g. a reference not yet covered) resolve to
    // nothing and PANTHER_API falls back to its prior, unfiltered behavior.
    def panther_prefix_map = new groovy.json.JsonSlurper()
                    .parse(new File("${projectDir}/assets/panther_species_prefixes.json"))
    ch_ref_panther_prefixes = ch_refs
                    .map { meta, gtf, fasta -> meta.id.toString() }
                    .collect()
                    .map { taxa -> taxa.collect { panther_prefix_map[it]?.prefix }.findAll { it }.unique().join(',') }

    ch_input_fa  = Channel.of(params.input)
                        .map { [[id: params.project_id], it] }
    ch_input_gtf = Channel.of(params.gtf)
                        .map { [[id: params.project_id], it] }
    ch_input     = ch_input_gtf.join(ch_input_fa)

    GUNZIP_REF( ch_refs )
    GUNZIP_IN( ch_input )


    GFFREAD(
        GUNZIP_REF.out.gtf.map { it[1] },
        GUNZIP_REF.out.fasta)

    ch_fasta = Channel.empty()

    GFFREAD_1(
        GUNZIP_IN.out.gtf.map { it[1] },
        GUNZIP_IN.out.fasta)

    GFFREAD_1.out.transcripts.join(GUNZIP_IN.out.gtf)
        .set { ch_query_tx_gtf }

    ch_query_tx_gtf.view()

    TRANSDECODER(
        ch_query_tx_gtf,
        params.pfam,
        params.project_id
    )
    .peptide_fasta.set { ch_fasta }


    ch_ref_tx  = GFFREAD.out.transcripts
    ch_ref_gtf = GUNZIP_REF.out.gtf
    ch_ref_tx.join(ch_ref_gtf)
            .set { ch_ref_tx_gtf }

    REF_TRANSDECODER(
        ch_ref_tx_gtf,
        params.pfam,
        params.project_id
    )

    ch_formatted_refs = Channel.empty()


    PRE_PROC(
        REF_TRANSDECODER.out.peptide_fasta.map { [it[0], it[1]] },
        REF_TRANSDECODER.out.peptide_fasta.map { [it[0], it[2]] }
    )


    PREP_INPUT(
        PRE_PROC.out,
        ch_fasta.first().map { [it[0], it[1]] },
        params.project_id
    )

    ch_ortho_f   = ORTHOFINDER(
                        PREP_INPUT.out.ortho_f,
                        params.project_id
                    )
    ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

    ch_diamond = Channel.empty()
    if (params.diamond_custom){
        DIAMOND(
            REF_TRANSDECODER.out.peptide_fasta.map { [it[0], it[1]] }
        )
        ch_diamond = DIAMOND.out.egg
    }

    ch_egg       = EGGNOGMAPPER(
                        PREP_INPUT.out.egg.first(),
                        params.project_id,
                        ch_diamond.ifEmpty([])
                    )
    ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions)


    ch_tree      = TREEGRAFTER(
                        PREP_INPUT.out.tree.first(),
                        params.sup_data_panther,
                        params.project_id
                    )
    //ch_versions.mix(TREEGRAFTER.out.versions)
    ch_panther   = PANTHER_API(ch_tree.splitText(by: 10, file: 'chunk.out'), ch_ref_panther_prefixes)
    COLLECT_CHUNKS(ch_panther.collect())


    ch_ortho_l   = ORTHOLOGER(
                        PREP_INPUT.out.ortho_l_data,
                        PREP_INPUT.out.ortho_l_work,
                        PREP_INPUT.out.ortho_l_work.first().map { it[1] },
                        params.project_id
                    )
    ch_versions = ch_versions.mix(ORTHOLOGER.out.versions)



    // -------------------------------------------------------------------------
    // SYNTENY (i-ADHoRe) — runs once per reference species; feeds into
    // consensus voting as the 5th evidence source.  Disable with --run_synteny false.
    // -------------------------------------------------------------------------

    ch_synteny = Channel.empty()
    if (params.run_synteny) {
        // Build per-scaffold gene-order list files for query and each reference
        GTF_TO_GENELISTS_QUERY(
            TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, gtf ] }
        )
        GTF_TO_GENELISTS_REF(
            REF_TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, gtf ] }
        )

        // DIAMOND all-vs-all blast table: query proteins vs each reference (one job per ref)
        ch_query_pep = TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, pep ] }
        ch_ref_pep   = REF_TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, pep ] }
        DIAMOND_BLASTP(ch_query_pep.first(), ch_ref_pep)

        // Join reference genelists with their blast table (both keyed by ref_meta)
        ch_ref_iadhore = GTF_TO_GENELISTS_REF.out.genelists
            .join(DIAMOND_BLASTP.out.blast)
        // shape: [ ref_meta, genelists_dir, blast_table ]

        // Broadcast single query genelists against every reference
        IADHORE(
            GTF_TO_GENELISTS_QUERY.out.genelists.first(),
            ch_ref_iadhore
        )
        MULTIPLICON_TO_GENES(IADHORE.out.outdir)
        SYNTENY_POST(IADHORE.out.synteny, params.sup_ensembl_data)
        ch_synteny  = SYNTENY_POST.out.synteny
        ch_versions = ch_versions.mix(IADHORE.out.versions)
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions.first())
        ch_versions = ch_versions.mix(GTF_TO_GENELISTS_QUERY.out.versions.first())
    }

    // -------------------------------------------------------------------------
    // Structural comparison (foldseek). Default evidence is ProstT5's
    // predicted 3Di sequence (no folding at all); fold_method='alphafold'
    // switches to user-supplied PDB structures in query_pdb_dir/ref_pdb_dir,
    // falling back to ProstT5 for any protein missing one. Standalone
    // diagnostic output -- not wired into CONSENSUS_VOTE.
    if (params.run_structural) {
        if (!params.prostt5_weights_dir) {
            error "params.prostt5_weights_dir must be set when params.run_structural is true -- a persistent (S3/JuiceFS) location for the one-time ProstT5 weight download, reused across runs instead of re-downloaded."
        }

        // Skip the (multi-GB) ProstT5 weight download if a previous run
        // already populated this persistent location.
        ch_weights_marker = file("${params.prostt5_weights_dir}/.download_complete")
        if (ch_weights_marker.exists()) {
            ch_prostt5_weights = Channel.value(file("${params.prostt5_weights_dir}/weights"))
        } else {
            FOLDSEEK_DOWNLOAD_PROSTT5()
            ch_prostt5_weights = FOLDSEEK_DOWNLOAD_PROSTT5.out.weights
            ch_versions = ch_versions.mix(FOLDSEEK_DOWNLOAD_PROSTT5.out.versions)
        }

        ch_query_pep_struct = TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, pep ] }
        ch_ref_pep_struct    = REF_TRANSDECODER.out.peptide_fasta.map { meta, pep, gtf -> [ meta, pep ] }

        // A real directory is only passed through in alphafold mode; the
        // bundled empty assets/NO_FILE placeholder means "not supplied" --
        // FOLDSEEK_CREATEDB tells the two apart with `[ -d ... ]`, not by name.
        ch_query_pdb_dir = (params.fold_method == 'alphafold' && params.query_pdb_dir)
            ? Channel.value(file(params.query_pdb_dir))
            : Channel.value(file("${projectDir}/assets/NO_FILE"))
        ch_ref_pdb_dir = (params.fold_method == 'alphafold' && params.ref_pdb_dir)
            ? Channel.value(file(params.ref_pdb_dir))
            : Channel.value(file("${projectDir}/assets/NO_FILE"))

        FOLDSEEK_CREATEDB_QUERY(
            ch_query_pep_struct,
            ch_prostt5_weights.first(),
            ch_query_pdb_dir.first(),
            params.fold_method
        )
        FOLDSEEK_CREATEDB_REF(
            ch_ref_pep_struct,
            ch_prostt5_weights.first(),
            ch_ref_pdb_dir.first(),
            params.fold_method
        )
        FOLDSEEK_SEARCH(
            FOLDSEEK_CREATEDB_QUERY.out.db.first(),
            FOLDSEEK_CREATEDB_REF.out.db
        )

        ch_versions = ch_versions.mix(FOLDSEEK_CREATEDB_QUERY.out.versions.first())
        ch_versions = ch_versions.mix(FOLDSEEK_CREATEDB_REF.out.versions.first())
        ch_versions = ch_versions.mix(FOLDSEEK_SEARCH.out.versions.first())
    }

    // -------------------------------------------------------------------------

    ch_ortho_f   = ch_ortho_f.ortho_f.ifEmpty(PREP_INPUT.out.blank)
                                    .branch {
                                        ortho_f: it
                                        na     : !it
                                    }

    ch_out_dir   = STAGE_OUTS(
                        ch_ortho_f.ortho_f,
                        ch_ortho_l.loger,
                        ch_egg.egg.first(),
                        COLLECT_CHUNKS.out.first()
                    )


    POST_PROC(
        ch_out_dir,
        PREP_INPUT.out.ortho_l_data,
        PREP_INPUT.out.egg.first(),
        params.sup_ensembl_data,
        params.taxa_db
    )


    CONSENSUS_VOTE(
        POST_PROC.out.ortho_l,
        POST_PROC.out.ortho_f,
        POST_PROC.out.eggnog.first(),
        POST_PROC.out.tree.first(),
        ch_synteny,
        ch_fasta.first().map { it[1] }
    )


    /* CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml')) */

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
