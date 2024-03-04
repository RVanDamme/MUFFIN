process chopper {
    label 'chopper'
    errorStrategy = { task.exitStatus in [14, -1] ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/quality_control/nanopore/", mode: 'copy', pattern: "*_cleaned.fastq" 
    input:
    tuple val(name), path(part)
    output:
    tuple val(name), path("${name}_cleaned.fastq")
    shell:
    """
    cat !{part} | chopper -q 10 -l ${params.short_qc} -t ${task.cpus} --headcrop 100 > ${name}_cleaned.fastq
    """
}

// process merge {
//     label 'ubuntu'
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     publishDir "${params.output}/${name}/assemble/quality_control/nanopore/", mode: 'copy', pattern: "*_all.fastq" 
//     input:
//     tuple val(name) , path(filtered)
//     output:
//     tuple val(name), path("${name}_all.fastq")
//     script:
//     """
//     cat *.fastq > ${name}_all.fastq
//     """

// }
