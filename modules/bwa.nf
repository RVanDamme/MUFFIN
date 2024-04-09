process bwa {
    label 'bwa'

    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.9"
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "illumina.bam"  
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(assembly), path(illumina)
    output:
    tuple val(name) , path("illumina_sorted.bam")
    script:
    """
    bwa index -p illumina -a bwtsw ${assembly}
    bwa mem illumina ${illumina[0]} ${illumina[1]} -t ${task.cpus} > illumina.sam
    samtools view -bS illumina.sam > illumina.bam
    samtools sort -@ ${task.cpus} -o illumina_sorted.bam illumina.bam
    rm illumina.*
    """
}

process extra_bwa {
    label 'bwa'

    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.9"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "illumina.bam"  
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(illumina)
    output:
    tuple val(name) , path("*_sorted.bam")
    script:
    """
    bwa index -p illumina -a bwtsw ${assembly}
    bwa mem illumina ${illumina} -t ${task.cpus} > illumina.sam
    samtools view -bS illumina.sam > illumina.bam
    samtools sort -@ ${task.cpus} -o ${illumina[0]}_sorted.bam illumina.bam
    rm illumina.*
    """
}

process bwa_bin {
    label 'bwa'

    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.9"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    //publishDir "${params.output}/${name}_bam/", mode: 'copy', pattern: "illumina.bam"  
    //SINCE THIS module is use multiple times it migh not be advise to output the same name file mutiple times
    input:
    tuple val(name), path(assembly), path(illumina)
    output:
    tuple val(name) , path("*_sorted.bam")
    script:
    """
    bin_id=\$(basename ${assembly} | sed -r "s/\\.\\w+//2")

    bwa index -p illumina -a bwtsw ${assembly}
    bwa mem illumina ${illumina[0]} ${illumina[1]} -t ${task.cpus} > illumina.sam
    samtools view -bS illumina.sam > illumina.bam
    samtools sort -@ ${task.cpus} -o ./bin_map/illumina/\$bin_id_sorted.bam illumina.bam
    rm illumina.*
    """
}