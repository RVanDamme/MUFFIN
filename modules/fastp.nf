process fastp {
    label 'fastp'
    publishDir "${params.output}/${name}/assemble/quality_control/illumina/", mode: 'copy', pattern: "*_R*_clean.fastq"
    input:
    tuple val(name), file(illumina)
    output:
    tuple val(name), file("*_R?_clean.fastq")
    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """
}

process fastp_rna {
    label 'fastp'
    publishDir "${params.output}/${name}/annotate/rna_quality_control/", mode: 'copy', pattern: "*_R*_clean.fastq"
    input:
    tuple val(name), file(illumina)
    output:
    tuple val(name), file("*_R?_clean.fastq")
    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """
}
