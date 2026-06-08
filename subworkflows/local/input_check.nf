//
// Validate the reference samplesheet and emit [ meta, gtf, fasta ] channels.
// Expected columns: taxa_id, gtf_file, fasta_file
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    Channel.fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { create_ref_channel(it) }
        .set { refs }

    emit:
    refs  // channel: [ val(meta), path(gtf), path(fasta) ]
}

def create_ref_channel(LinkedHashMap row) {
    if (!row.taxa_id) {
        error "ERROR: samplesheet -> 'taxa_id' column is missing or empty."
    }
    if (!file(row.gtf_file).exists()) {
        error "ERROR: samplesheet -> GTF file does not exist!\n${row.gtf_file}"
    }
    if (!file(row.fasta_file).exists()) {
        error "ERROR: samplesheet -> FASTA file does not exist!\n${row.fasta_file}"
    }

    return [ [ id: row.taxa_id ], file(row.gtf_file), file(row.fasta_file) ]
}
