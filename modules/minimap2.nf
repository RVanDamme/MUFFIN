process minimap2 {
    label 'minimap2'
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    set val(name), file(assembly), file(ont)
    output:
    set val(name) , file("ont_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -@ ${task.cpus} -o ont_sorted.bam ont.bam
    """
}

process minimap_polish {
    label 'minimap2'
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    set val(name), file(assembly), file(ont)
    output:
    set val(name) , file("ont.sam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    """
}

process extra_minimap2 {
    label 'minimap2'
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    set val(name), file(assembly), file(ont)
    output:
    set val(name) , file("*_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -@ ${task.cpus} -o ${ont}_sorted.bam ont.bam
    """
}


process minimap2_bin {
    label 'minimap2'
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
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