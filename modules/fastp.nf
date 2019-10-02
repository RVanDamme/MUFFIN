process fastp {
    label 'fastp'
    if (params.out_qc == true ) { publishDir "${params.output}/${name}_illumina_qc/", mode: 'copy', pattern: "${name}_R?_clean.fastq"}
    input:
    set val(name), file(illumina)

    output:
    set val(name), file(["${name}_R1_clean.fastq","${name}_R2_clean.fastq"])

    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """

}