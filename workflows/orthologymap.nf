
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

//Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

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
include { BLASTP } from "../modules/local/blastp.nf"
include { DAGCHAINER } from "../modules/local/dagchainer.nf"

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
    
    ch_input_fa  = Channel.of(params.input).map { [[id: params.project_id], it] } 
    ch_input_gtf = Channel.of(params.gtf).map { [[id: params.project_id], it] }
    ch_input     = ch_input_gtf.join(ch_input_fa)
    
    GUNZIP_REF( ch_refs )
    GUNZIP_IN( ch_input )

        
    GFFREAD(
        GUNZIP_REF.out.gtf.map { it[1] }, 
        GUNZIP_REF.out.fasta)

    ch_fasta = Channel.empty() 
                 
    ch_query_tx  = GUNZIP_IN.out.fasta
    ch_query_gtf = GUNZIP_IN.out.gtf
    ch_query_tx.join(ch_query_gtf).set { ch_query_tx_gtf } 
    
    
    TRANSDECODER(
        ch_query_tx_gtf,
        params.pfam, 
        params.project_id
    )
    .peptide_fasta.set { ch_fasta }

    
    ch_ref_tx  = GFFREAD.out.transcripts
    ch_ref_gtf = GUNZIP_REF.out.gtf
    ch_ref_tx.join(ch_ref_gtf).set { ch_ref_tx_gtf } 

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
    ch_versions.mix(ORTHOFINDER.out.versions)


    ch_egg       = EGGNOGMAPPER(
                       PREP_INPUT.out.egg.first(), 
                       params.project_id
                   )
    ch_versions.mix(EGGNOGMAPPER.out.versions)

    
    ch_tree      = TREEGRAFTER(
                       PREP_INPUT.out.tree.first(), 
                       params.sup_data_panther, 
                       params.project_id
                   )
    //ch_versions.mix(TREEGRAFTER.out.versions)
    ch_panther   = PANTHER_API(ch_tree.splitText(by: 10, file: 'chunk.out'))
    COLLECT_CHUNKS(ch_panther.collect())

    
    ch_ortho_l   = ORTHOLOGER(
                       PREP_INPUT.out.ortho_l_data, 
                       PREP_INPUT.out.ortho_l_work, 
                       PREP_INPUT.out.ortho_l_work.first().map { it[1] }, 
                       params.project_id
                   )
    ch_versions.mix(ORTHOLOGER.out.versions)


    
    //ch_blastp  = BLASTP(
    //                 ch_fasta.first().map { [it[0], it[1]] },
    //                 PRE_PROC.out
    //             )
  
    //ch_dag     = DAGCHAINER(
    //                 ch_blastp.tbl,
    //                 ch_query_gtf.first(),
    //                 ch_ref_gtf, 
    //                 params.project_id
    //             )

    ch_ortho_f   = ch_ortho_f.ortho_f.ifEmpty(PREP_INPUT.out.blank).branch {
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
        POST_PROC.out.tree.first()
    )
                    
    
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

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
