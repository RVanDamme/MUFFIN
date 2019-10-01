process minimap2 {
    label 'minimap2'
    publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    input:
    set val(name), file(assembly), file(ont)
    output:
    set val(name) , file("ont_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -o ont_sorted.bam ont.bam
    """
}