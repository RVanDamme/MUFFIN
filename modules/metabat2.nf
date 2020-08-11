process metabat2 {
    maxForks 1
    label 'metabat2'
    publishDir "${params.output}/${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(assembly), path(ont_bam), path(illumina_bam)
    output:
    tuple val(name), path("bins_dir")
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
    metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bins -t ${task.cpus}
    """
}

process metabat2_extra {
    maxForks 1
    label 'metabat2'
    publishDir "${params.output}/${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir" 
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(assembly), path(ont_bam), path(illumina_bam)
    path(extra_bam)
    output:
    tuple val(name), path("bins_dir")
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
    metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bin -t ${task.cpus}
    """

}