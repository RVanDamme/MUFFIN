process discard_short {
    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name) , path(part)
    output:
    tuple val(name), path("filtered_${part}")
    shell:
    """
        cat !{part} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${params.short_qc}' | sed 's/\\t/\\n/g' > "filtered_${part}"

    """
}



process filtlong {
    label 'filtlong'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name) , path(filtered)
    output:
    tuple val(name) , path("clean_${filtered}")
    script:
    """
    filtlong --min_length ${params.short_qc} --keep_percent 90 --target_bases 500000000 ${filtered} > clean_${filtered}
    """
}


process merge {
    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/quality_control/nanopore/", mode: 'copy', pattern: "*_all.fastq" 
    input:
    tuple val(name) , path(filtered)
    output:
    tuple val(name), path("${name}_all.fastq")
    script:
    """
    cat *.fastq > ${name}_all.fastq
    """

}