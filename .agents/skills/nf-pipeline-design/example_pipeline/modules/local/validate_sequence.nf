nextflow.enable.types = true

process VALIDATE_SEQUENCE {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    tuple(meta: Map, fasta: Path)

    output:
    validated = tuple(meta, file("*.validated.fa"))

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    # Strip any non-ACGT characters from the sequence body so downstream
    # composition maths works on a clean alphabet. Headers are preserved.
    awk 'BEGIN { IGNORECASE=1 }
         /^>/ { print; next }
         { gsub(/[^ACGTacgt]/, ""); print toupper(\$0) }' ${fasta} > ${prefix}.validated.fa
    """
}
