//
// SPLIT_INPUT_FASTA: split the input FASTA into per-record files and attach a sample
// meta map to each. Emits one channel item per sequence record.
//
include { SPLIT_FASTA } from '../../modules/local/split_fasta'

workflow SPLIT_INPUT_FASTA {
    take:
    ch_fasta

    main:
    ch_split = SPLIT_FASTA(ch_fasta)

    ch_records = ch_split
        .flatten()
        .map { record -> tuple([id: record.baseName], record) }

    emit:
    records = ch_records
}
