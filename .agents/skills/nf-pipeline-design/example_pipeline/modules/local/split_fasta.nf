nextflow.enable.types = true

process SPLIT_FASTA {
    tag "${fasta.baseName}"
    label 'process_single'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    fasta: Path

    output:
    records: Set<Path> = files("records/*.fa")

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    """
    mkdir -p records
    awk 'BEGIN { out=""; seq="" }
         /^>/ {
             if (out) { printf "%s", seq > out; close(out); seq="" }
             name=\$1
             gsub(/^>/, "", name)
             out="records/" name ".fa"
             print \$0 > out
             next
         }
         { seq = seq \$0 "\\n" }
         END { if (out) { printf "%s", seq > out; close(out) } }' ${fasta}
    """
}
