nextflow.enable.types = true

process COMPUTE_GC_CONTENT {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    tuple(meta: Map, fasta: Path)

    output:
    gc = tuple(meta, file("*.gc.tsv"))

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk 'BEGIN { gc=0; total=0 }
         /^>/ { next }
         { for (i=1; i<=length(\$0); i++) {
               c = toupper(substr(\$0, i, 1))
               if (c == "G" || c == "C") gc++
               if (c == "A" || c == "C" || c == "G" || c == "T") total++
             }
         }
         END {
             if (total == 0) { ratio = 0 } else { ratio = gc/total }
             printf "sample\\tgc\\ttotal\\tratio\\n${meta.id}\\t%d\\t%d\\t%.4f\\n", gc, total, ratio
         }' ${fasta} > ${prefix}.gc.tsv
    """
}
