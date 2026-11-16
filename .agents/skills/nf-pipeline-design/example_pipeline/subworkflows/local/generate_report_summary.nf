//
// GENERATE_REPORT_SUMMARY: collate per-sequence GC and k-mer metric TSVs into a
// single summary TSV. Input-shaping (dropping meta and collecting across sequences)
// lives here, so the module stays a simple "take TSVs, emit summary" tool.
//
include { GENERATE_SUMMARY } from '../../modules/local/generate_summary'

workflow GENERATE_REPORT_SUMMARY {
    take:
    ch_gc        // tuple(meta, gc_tsv)
    ch_kmer      // tuple(meta, kmer_tsv)

    main:
    ch_gc_files   = ch_gc  .map { _meta, tsv -> tsv }.collect()
    ch_kmer_files = ch_kmer.map { _meta, tsv -> tsv }.collect()

    summary_result = GENERATE_SUMMARY(ch_gc_files, ch_kmer_files)

    emit:
    summary = summary_result.summary
}
