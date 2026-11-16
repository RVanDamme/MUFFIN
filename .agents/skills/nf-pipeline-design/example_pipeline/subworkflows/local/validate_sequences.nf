//
// VALIDATE_SEQUENCES: filter each per-sequence record down to a clean ACGT alphabet.
//
include { VALIDATE_SEQUENCE } from '../../modules/local/validate_sequence'

workflow VALIDATE_SEQUENCES {
    take:
    ch_records     // tuple(meta, fasta_record)

    main:
    validate_result = VALIDATE_SEQUENCE(ch_records)

    emit:
    sequences = validate_result.validated
}
