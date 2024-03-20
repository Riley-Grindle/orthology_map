
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


include { TRANSDECODER } from "../subworkflows/local/transdecoder.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ORTHOLOGYMAP {
    
    ch_versions = Channel.empty()

    ch_fasta = TRANSDECODER(params.input, params.gtf, params.pfam, params.project_id).peptide_fasta
    ch_ref_fastas = Channel.fromPath(params.aws_reference_ortho)

    ch_formatted_refs = PRE_PROC( ch_ref_fastas )

    PREP_INPUT(ch_formatted_refs, ch_fasta,  params.project_id)

    ch_ortho_f = ORTHOFINDER(PREP_INPUT.out.ortho_f.first(), params.project_id)
    ch_versions.mix(ORTHOFINDER.out.versions)    

    ch_egg = EGGNOGMAPPER(PREP_INPUT.out.egg.first(), params.project_id)
    ch_versions.mix(EGGNOGMAPPER.out.versions)

    ch_tree = TREEGRAFTER(PREP_INPUT.out.tree.first(), params.sup_data_panther, params.project_id)
   
    ch_ortho_l = ORTHOLOGER(PREP_INPUT.out.ortho_l_data.first(), PREP_INPUT.out.ortho_l_work.first(), PREP_INPUT.out.ortho_l_work.first(), params.project_id)
    ch_versions.mix(ORTHOLOGER.out.versions)

    ch_panther = PANTHER_API(ch_tree.splitText(by: 10, file: 'chunk.out'))
    COLLECT_CHUNKS(ch_panther.collect())
  
    ch_ortho_f = ch_ortho_f.ortho_f.ifEmpty(PREP_INPUT.out.blank).branch { 
									ortho_f: it    
									na     : !it
								 }
    
    ch_out_dir = STAGE_OUTS(ch_ortho_f.ortho_f, ch_ortho_l.loger, ch_egg.egg, COLLECT_CHUNKS.out.first())

    ch_outs = POST_PROC(ch_out_dir, PREP_INPUT.out.ortho_l_data.first(), PREP_INPUT.out.egg.first(), params.sup_ensembl_data)
    
    CUSTOM_DUMPSOFTWAREVERSIONS (ch_versions.unique().collectFile(name: 'collated_versions.yml'))
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
