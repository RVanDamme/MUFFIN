nextflow.enable.types = true

process COMPUTE_KMER_PROFILE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    tuple(meta: Map, fasta: Path)
    kmer_size: Integer

    output:
    kmer = tuple(meta, file("*.kmer.tsv"))

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk -v k=${kmer_size} 'BEGIN { IGNORECASE=1 }
        /^>/ { next }
        { seq = seq toupper(\$0) }
        END {
            for (i = 1; i + k - 1 <= length(seq); i++) {
                kmer = substr(seq, i, k)
                counts[kmer]++
            }
            printf "kmer\\tcount\\n"
            for (kmer in counts) printf "%s\\t%d\\n", kmer, counts[kmer]
        }' ${fasta} | sort > ${prefix}.kmer.tsv
    """
}
