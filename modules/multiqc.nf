process chopper {
    label 'multiqc'

    errorStrategy = { task.exitStatus in [14, -1] ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/report/", mode: 'copy', pattern: "*" 
    input:
    tuple val(name), path(ont)
    output:
    tuple val(name), path("${name}_cleaned.fastq")
    shell:
    """
    chopper -q 10 -l ${params.short_qc} -t ${task.cpus} --headcrop 100 -i ${ont} > ${name}_cleaned.fastq
    """
}