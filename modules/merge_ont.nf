process merge_ont {
    label 'ubuntu'
    publishDir "${params.output}/${name}_clean/", mode: 'copy', pattern: "*_clean.fastq"
    input:
    set val(name1) , file(part1)
    set val(name2) , file(part2)
    set val(name3) , file(part3)
    set val(name4) , file(part4)
    output:
    set val(name1), file("${name1}_clean.fastq")
    script:
    """"
    cat ${part1} ${part2} ${part3} ${part4} > ${name1}_clean.fastq
    """"

}