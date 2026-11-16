//
// ANALYZE_GC: compute the GC fraction of each sequence record.
//
include { COMPUTE_GC_CONTENT } from '../../modules/local/compute_gc_content'

workflow ANALYZE_GC {
    take:
    ch_sequences   // tuple(meta, fasta)

    main:
    gc_result = COMPUTE_GC_CONTENT(ch_sequences)

    emit:
    gc = gc_result.gc
}
