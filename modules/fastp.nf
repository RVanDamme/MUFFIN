process fastp {
    label 'fastp'

    conda 'bioconda::fastp=0.20.0'

    publishDir "${params.output}/${name}/assemble/quality_control/illumina/", mode: 'copy', pattern: "*_R*_clean.fastq"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(illumina)
    output:
    tuple val(name), path("*_R?_clean.fastq")
    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """
}

process fastp_rna {
    label 'fastp'

    conda 'bioconda::fastp=0.20.0'

    publishDir "${params.output}/${name}/annotate/rna_quality_control/", mode: 'copy', pattern: "*_R*_clean.fastq"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(illumina)
    output:
    tuple val(name), path("*_R?_clean.fastq")
    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """
}
