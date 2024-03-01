// process metabat2 {
//     maxForks 1
//     label 'metabat2'
//     publishDir "${params.output}/${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     tuple val(name), path(assembly), path(ont_bam), path(illumina_bam)
//     output:
//     tuple val(name), path("bins_dir")
//     script:
//     """
//     jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
//     metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bins -t ${task.cpus}
//     """
// }

// process metabat2_extra {
//     maxForks 1
//     label 'metabat2'
//     publishDir "${params.output}/${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir" 
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     tuple val(name), path(assembly), path(ont_bam), path(illumina_bam)
//     path(extra_bam)
//     output:
//     tuple val(name), path("bins_dir")
//     script:
//     """
//     jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
//     metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bin -t ${task.cpus}
//     """

// }

// process metabat2 {
//     maxForks 1
//     label 'metabat2'
//     publishDir "\${params.output}/\${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
    
//     output:
//     tuple val(name), path("bins_dir")
    
//     script:
//     """
//     # Initialize an empty file (concatenate depths)
//     > depth.txt

//     # Check if ont_bam is provided and process it
//     if [ ! -z "${ont_bam}" ]; then
//         jgi_summarize_bam_contig_depths --outputDepth ont_depth.txt ${ont_bam}
//         cat ont_depth.txt >> depth.txt
//     fi

//     # Check if illumina_bam is provided and process it
//     if [ ! -z "${illumina_bam}" ]; then
//         jgi_summarize_bam_contig_depths --outputDepth illumina_depth.txt ${illumina_bam}
//         cat illumina_depth.txt >> depth.txt
//     fi

//     # Check if extra_bam is provided and process it
//     if [ ! -z "${extra_bam}" ]; then
//         jgi_summarize_bam_contig_depths --outputDepth illumina_depth.txt ${extra_bam}
//         cat illumina_depth.txt >> depth.txt
//     fi

//     # Run MetaBAT2 if depth.txt is not empty
//     if [ -s depth.txt ]; then
//         metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bins -t ${task.cpus}
//     else
//         echo "No BAM files provided for depth calculation."
//         exit 1
//     fi
//     """
// }


process metabat2 {
    maxForks 1
    label 'metabat2'
    publishDir "\${params.output}/\${name}/assemble/binning/metabat2/", mode: 'copy', pattern: "bins_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(assembly), path(bam)
    
    output:
    tuple val(name), path("bins_dir")
    
    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    metabat2 -i ${assembly} -a depth.txt -o bins_dir/metabat_bins -t ${task.cpus}
    """
}