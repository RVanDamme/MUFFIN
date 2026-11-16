nextflow.enable.types = true

process GENERATE_SUMMARY {
    tag "composition_summary"
    label 'process_single'

    container 'quay.io/nf-core/ubuntu:20.04'

    input:
    gc_files: List<Path>
    kmer_files: List<Path>

    stage:
    stageAs gc_files,   'gc/*'
    stageAs kmer_files, 'kmer/*'

    output:
    summary: Path = file("composition_summary.tsv")

    topic:
    tuple('awk', eval('awk --version | head -1')) >> 'versions'

    script:
    """
    # Stitch GC rows across samples (keep a single header).
    awk 'FNR==1 && NR!=1 { next } { print }' gc/*.tsv > composition_summary.tsv

    # Append richness: number of distinct k-mers observed per sample.
    printf "\\nsample\\tdistinct_kmers\\n" >> composition_summary.tsv
    for f in kmer/*.tsv; do
        sample=\$(basename "\$f" .kmer.tsv)
        distinct=\$(awk 'NR>1' "\$f" | wc -l | tr -d ' ')
        printf "%s\\t%s\\n" "\$sample" "\$distinct" >> composition_summary.tsv
    done
    """
}
