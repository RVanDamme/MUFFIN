process metabat2 {
    maxForks 1
    label 'metabat2'
    publishDir "${params.output}/${name}/metabat2_bins/", mode: 'copy', pattern: "bins_dir"
    input:
    tuple val(name), file(assembly), file(ont_bam), file(illumina_bam)
    output:
    tuple val(name), file("bins_dir")
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
    metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bins -t ${task.cpus}
    """
}

process metabat2_extra {
    maxForks 1
    label 'metabat2'
    publishDir "${params.output}/${name}/metabat2_bins/", mode: 'copy', pattern: "bins_dir" 
    input:
    tuple val(name), file(assembly), file(ont_bam), file(illumina_bam)
    file(extra_bam)
    output:
    tuple val(name), file("bins_dir")
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
    metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bin -t ${task.cpus}
    """

}