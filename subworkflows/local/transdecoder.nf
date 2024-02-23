
include { GTF_2_GENETX_MAP } from '../../modules/local/gtf_2_genetx_map.nf'
include { TRANSDECODER_LONGORF } from '../../modules/local/longorf.nf'
include { HMMER_HMMSCAN } from '../../modules/local/hmmscan.nf'
include { TRANSDECODER_PREDICT } from '../../modules/local/predict.nf'
include { GUNZIP as GUNZIP_FASTA } from '../../modules/local/gunzip/main'
include { GUNZIP as GUNZIP_GTF } from '../../modules/local/gunzip/main'
include { GUNZIP as GUNZIP_PFAM } from '../../modules/local/gunzip/main'


workflow TRANSDECODER {

    take:
    fasta                  // file: path/to/fasta
    gtf                    // file: path/to/gtf
    pfam                   // file: path/to/pfam_hmm
    project_id             // val : "project_id"

    main:
    
    ch_versions = Channel.empty()
    

    // Unzip fasta
   
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    // Unzip gtf 

    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
        } else {
            ch_gtf = Channel.value(file(gtf))
    }
    }
    
    // Unzip pfam hmm file
   
    if (pfam) {
	if (pfam.endsWith('.gz')) {
            ch_pfam	= GUNZIP_PFAM ( [ [:], pfam ] ).gunzip.map { it[1] }
        } else {
            ch_pfam = Channel.value(file(pfam))
    }
    }

    // Create gene to transcript map file to be used in LongOrfs
    
    GTF_2_GENETX_MAP(ch_gtf)

    //Run Transdecoder.LongOrf

    TRANSDECODER_LONGORF(ch_fasta, GTF_2_GENETX_MAP.out.genetx_map)
    
    // Build Domain hit table from longest orfs

    HMMER_HMMSCAN(ch_pfam, TRANSDECODER_LONGORF.out.pep)
    
    // Run Transdecoder predict while preserving protein domain hits from pfam db

    ch_pep = TRANSDECODER_PREDICT(ch_fasta, HMMER_HMMSCAN.out.table, TRANSDECODER_LONGORF.out.out_dir, project_id).pep 
    

     emit:
     peptide_fasta = ch_pep 
     versions      = ch_versions   



}

