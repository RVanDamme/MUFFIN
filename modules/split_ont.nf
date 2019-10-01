process split_ont {
    label 'ubuntu'
    input:
        set val(name), file(ont)
    output:
    val(name), file("part0.fastq")    emit:   ont_1
    val(name), file("part1.fastq")    emit:   ont_2
    val(name), file("part2.fastq")    emit:   ont_3
    val(name), file("part3.fastq")    emit:   ont_4
    script:
    """
    total_line=$(wc -l !{ont} | cut -d ' ' -f 1)
    line_per_file=$((\$total_line/4))
    split -d -a 1 -l \$line_per_file !{ont} part
    for part in part*; do mv \$part \$part".fastq"; done
    """

}