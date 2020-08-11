process minimap2 {
    label 'minimap2'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(ont)
    output:
    tuple val(name) , path("ont_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -@ ${task.cpus} -o ont_sorted.bam ont.bam
    rm ont.*
    """
}

process minimap_polish {
    label 'minimap2'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(ont)
    output:
    tuple val(name) , path("ont.paf")
    script:
    """
    minimap2 -x map-ont ${assembly} ${ont} > ont.paf
    """
}

process extra_minimap2 {
    label 'minimap2'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(ont)
    output:
    tuple val(name) , path("*_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -@ ${task.cpus} -o ${ont}_sorted.bam ont.bam
    rm ont.*
    """
}


process minimap2_bin {
    label 'minimap2'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "ont.bam"
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(ont)
    output:
    tuple val(name) , path("ont_sorted.bam")
    script:
    """
    minimap2 -ax map-ont ${assembly} ${ont} > ont.sam
    samtools view -bS ont.sam > ont.bam
    samtools sort -o ont_sorted.bam ont.bam
    rm ont.*
    """
}