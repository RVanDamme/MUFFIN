process bam_merger {
    label 'samtools_merger'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5

    input:
    tuple val(name), path(bam1), path(bam2)

    output:
    tuple val(name) , path("merged_bam.bam")

    script:
    """
    samtools merge -@ \${task.cpus} -o \${name}_merged.bam \${bam1} \${bam2} 
    """
}