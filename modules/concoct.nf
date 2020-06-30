process concoct {
    maxForks 1
    label 'concoct'
    publishDir "${params.output}/${name}/concoct_bins/", mode: 'copy', pattern: "fasta_bins"
    input:
    tuple val(name), file(assembly), file(ont_bam), file(illumina_bam)
    output:
    tuple val(name), file("fasta_bins")
    script:
    """
    mkdir concoct_out
    cut_up_fasta.py ${assembly} -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
    samtools index -@ ${task.cpus} ${ont_bam}
    samtools index -@ ${task.cpus} ${illumina_bam}
    concoct_coverage_table.py contigs_10K.bed *.bam > coverage_table.tsv
    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_out --thread ${task.cpus}
    merge_cutup_clustering.py concoct_out/clustering_gt1000.csv > concoct_out/clustering_merged.csv
    mkdir fasta_bins
    extract_fasta_bins.py ${assembly} concoct_out/clustering_merged.csv --output_path fasta_bins

    """

}

process concoct_extra {
    maxForks 1
    label 'concoct'
    publishDir "${params.output}/${name}/assemble/binning/concoct/", mode: 'copy', pattern: "fasta_bins"
    input:
    tuple val(name), path(assembly), path(ont_bam), path(illumina_bam)
    path(extra_bam)
    output:
    tuple val(name), path("fasta_bins")
    script:
    """
    mkdir concoct_out
    cut_up_fasta.py ${assembly} -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
    samtools index -@ ${task.cpus} ${ont_bam}
    samtools index -@ ${task.cpus} ${illumina_bam}
    ls ${extra_bam} | xargs -n1 -P5 samtools index  -@ ${task.cpus}
    concoct_coverage_table.py contigs_10K.bed *.bam > coverage_table.tsv
    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_out --thread ${task.cpus}
    merge_cutup_clustering.py concoct_out/clustering_gt1000.csv > concoct_out/clustering_merged.csv
    mkdir fasta_bins
    extract_fasta_bins.py ${assembly} concoct_out/clustering_merged.csv --output_path fasta_bins

    """

}