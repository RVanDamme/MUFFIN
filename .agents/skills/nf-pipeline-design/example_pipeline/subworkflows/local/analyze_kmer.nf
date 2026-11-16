//
// ANALYZE_KMER: compute the k-mer frequency profile of each sequence record.
//
include { COMPUTE_KMER_PROFILE } from '../../modules/local/compute_kmer_profile'

workflow ANALYZE_KMER {
    take:
    ch_sequences   // tuple(meta, fasta)
    kmer_size      // val

    main:
    kmer_result = COMPUTE_KMER_PROFILE(ch_sequences, kmer_size)

    emit:
    kmer = kmer_result.kmer
}
