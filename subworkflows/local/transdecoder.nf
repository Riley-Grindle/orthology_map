
include { GTF_2_GENETX_MAP       } from '../../modules/local/gtf_2_genetx_map.nf'
include { TRANSDECODER_LONGORF   } from '../../modules/local/longorf.nf'
include { HMMER_HMMSCAN          } from '../../modules/local/hmmscan.nf'
include { TRANSDECODER_PREDICT   } from '../../modules/local/predict.nf'
include { GUNZIP as GUNZIP_PFAM  } from '../../modules/local/gunzip/main'


workflow TRANSDECODER {

    take:
    tx_gtf      // ch  : tuple val(meta), path(.fa), path(.gtf)
    pfam        // file: path/to/pfam_hmm  (or null)
    project_id  // val : "project_id"

    main:

    ch_versions = Channel.empty()

    tx_gtf.map { meta, fasta, gtf -> [ meta, fasta ] }.set { ch_fasta }
    tx_gtf.map { meta, fasta, gtf -> [ meta, gtf   ] }.set { ch_gtf   }

    GTF_2_GENETX_MAP(ch_gtf)

    ch_fasta
        .join(GTF_2_GENETX_MAP.out.genetx_map)
        .multiMap { meta, fasta, genetx_map ->
            fasta_ch: [ meta, fasta ]
            map_ch:   [ meta, genetx_map ]
        }
        .set { ch_longorf_in }

    TRANSDECODER_LONGORF(ch_longorf_in.fasta_ch, ch_longorf_in.map_ch)
    ch_versions = ch_versions.mix(TRANSDECODER_LONGORF.out.versions)

    // Build domain-hit table only when a Pfam HMM database is supplied
    if (pfam) {
        if (pfam.endsWith('.gz')) {
            ch_pfam = GUNZIP_PFAM( [ [:], pfam ] ).gunzip.map { it[1] }
        } else {
            ch_pfam = Channel.value(file(pfam))
        }
        HMMER_HMMSCAN(ch_pfam, TRANSDECODER_LONGORF.out.pep)
        ch_pfam_table = HMMER_HMMSCAN.out.table
    } else {
        // Pass an empty list so TRANSDECODER_PREDICT skips --retain_pfam_hits
        ch_pfam_table = TRANSDECODER_LONGORF.out.pep.map { meta, pep -> [ meta, [] ] }
    }

    ch_pep = TRANSDECODER_PREDICT(
        TRANSDECODER_LONGORF.out.transcripts,
        ch_pfam_table,
        TRANSDECODER_LONGORF.out.out_dir,
        project_id
    ).pep

    ch_pep.join(ch_gtf).set { ch_combined }

    emit:
    peptide_fasta = ch_combined
    versions      = ch_versions

}
