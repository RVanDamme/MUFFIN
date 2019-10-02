process concoct {
    label 'concoct'
    label 'ubuntu'
    if (params.out_concoct == true) {publishDir "${params.output}/${name}_concoct/", mode: 'copy', pattern: "bins_dir/concoct_bin/fasta_bins/*"}
    input:
    set val(name), file(assembly), file(ont_bam), file(illumina_bam)
    output:
    set val(name), file("bins_dir/concoct_out/fasta_bins")
    script:
    """
    cut_up_fasta.py ${assembly} -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
    concoct_coverage_table.py contigs_10K.bed ${ont_bam},${illumina_bam} > coverage_table.tsv
    concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b bins_dir/concoct_out
    merge_cutup_clustering.py bins_dir/concoct_out/clustering_gt1000.csv > bins_dir/concoct_out/clustering_merged.csv
    mkdir bins_dir/concoct_out/fasta_bins
    extract_fasta_bins.py ${assembly} bins_dir/concoct_out/clustering_merged.csv --output_path bins_dir/concoct_bin/fasta_bins
    """

}