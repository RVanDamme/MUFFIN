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
    set val(name) , file("filtlong_${filtered}")
    script:
    """
    filtlong --min_length ${params.short_qc} --keep_percent 90 --target_bases 500000000 ${filtered} > clean_${filtered}
    """
}


process merge {
    label 'ubuntu'
    if (params.out_qc == true ) { publishDir "${params.output}/${name}_ont_qc/", mode: 'copy', pattern: "${name}.fastq" }
    input:
    set val(name) , file(filtered)
    output:
    set val(name), file("${name}.fastq")
    script:
    """
    cat *.fastq > ${name}.fastq
    """
}