process fastp {
    label 'fastp'
    if (params.out_qc == true ) { publishDir "${params.output}/${name}_illumina_qc/", mode: 'copy', pattern: "*_R?_clean.fastq"}
    input:
    set val(name), file(illumina)

    output:
    set val(name), file("*_R?_clean.fastq")

    script:
    """
    fastp -i ${illumina[0]} -I ${illumina[1]} -o ${name}_R1_clean.fastq -O ${name}_R2_clean.fastq
    """

}