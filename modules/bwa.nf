process bwa {
    label 'bwa'
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "illumina.bam"  
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    set val(name), file(assembly), file(illumina)
    output:
    set val(name) , file("illumina_sorted.bam")
    script:
    """
    bwa index -p illumina -a bwtsw ${assembly}
    bwa mem illumina ${illumina[0]} ${illumina[1]} -t ${task.cpus} > illumina.sam
    samtools view -bS illumina.sam > illumina.bam
    samtools sort -o illumina_sorted.bam illumina.bam
    """
}