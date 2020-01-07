process discard_short {
    label 'ubuntu'
    input:
    set val(name) , file(part)
    output:
    set val(name), file("filtered_${part}")
    shell:
    """
        cat !{part} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${params.short_qc}' | sed 's/\\t/\\n/g' > "filtered_${part}"

    """
}



process filtlong {
    label 'filtlong'
    input:
    set val(name) , file(filtered)
    output:
    set val(name) , file("clean_${filtered}")
    script:
    """
    filtlong --min_length ${params.short_qc} --keep_percent 90 --target_bases 500000000 ${filtered} > clean_${filtered}
    """
}


process merge {
    label 'ubuntu'
    publishDir "${params.output}/${name}/nanopore_qc_out/", mode: 'copy', pattern: "*_all.fastq" 
    input:
    set val(name) , file(filtered)
    output:
    set val(name), file("${name}_all.fastq")
    script:
    """
    cat *.fastq > ${name}_all.fastq
    """

}