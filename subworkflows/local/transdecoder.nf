
include { GTF_2_GENETX_MAP } from '../../modules/local/gtf_2_genetx_map.nf'
include { TRANSDECODER_LONGORF } from '../../modules/local/longorf.nf'
include { HMMER_HMMSCAN } from '../../modules/local/hmmscan.nf'
include { TRANSDECODER_PREDICT } from '../../modules/local/predict.nf'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/local/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from '../../modules/local/gunzip/main'
include { GUNZIP as GUNZIP_PFAM } from '../../modules/local/gunzip/main'


workflow TRANSDECODER {

    take:
    tx_gtf                 // ch  : tuple val(meta), path(.fa), path(.gtf)
    pfam                   // file: path/to/pfam_hmm
    project_id             // val : "project_id"

    main:
    
    ch_versions = Channel.empty()
    
    // Unzip pfam hmm file
   
    if (pfam) {
	if (pfam.endsWith('.gz')) {
            ch_pfam	= GUNZIP_PFAM ( [ [:], pfam ] ).gunzip.map { it[1] }
        } else {
            ch_pfam = Channel.value(file(pfam))
    }
    }

    // Create gene to transcript map file to be used in LongOrfs
    
    tx_gtf.map {
        meta, fasta, gtf ->
            return [meta, fasta]
    }
    .set { ch_fasta }
   
    tx_gtf.map {
        meta, fasta, gtf ->
            return [meta, gtf]
    }
    .set { ch_gtf }   
  
    GTF_2_GENETX_MAP(ch_gtf)

    //Run Transdecoder.LongOrf

    TRANSDECODER_LONGORF(ch_fasta, GTF_2_GENETX_MAP.out.genetx_map)
    
    // Build Domain hit table from longest orfs

    HMMER_HMMSCAN(ch_pfam, TRANSDECODER_LONGORF.out.pep)
    
    // Run Transdecoder predict while preserving protein domain hits from pfam db

    ch_pep = TRANSDECODER_PREDICT(TRANSDECODER_LONGORF.out.transcripts, HMMER_HMMSCAN.out.table, TRANSDECODER_LONGORF.out.out_dir, project_id).pep 
    
    ch_pep.join(ch_gtf).set { ch_combined }

     emit:
     peptide_fasta = ch_combined 
     versions      = ch_versions   



}

