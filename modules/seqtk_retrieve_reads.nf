process unmapped_retrieve {
    label 'seqtk'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/", mode: 'copy', pattern: "*unmapped_*.fastq"
    input:
    tuple val(name), path(bam), path(reads), val(sequencing_type)
    output:
    path("unmapped_*.fastq") optionnal true
    shell:
    """
    if [[ "${sequencing_type}" == "illumina" ]]; then
        ## illumina unmapped reads retrieval
        samtools view -f4 ${bam} > illumina_unmapped_contigs.sam
        cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
        seqtk subseq ${reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
        if [[ "${reads[1]}" != "" ]]; then
            seqtk subseq ${reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq
            fi
       
    elif [[ "${sequencing_type}" == "ont" ]]; then
        ## ONT unmapped reads retrieval
        samtools view -f4 ${bam} > ont_unmapped_contigs.sam
        cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
        seqtk subseq ${reads} ont_unmapped.list > unmapped_ONT.fastq
    fi

    rm *_unmapped_contigs.sam
    rm *_unmapped.list
    """

}

process unmapped_illumina_retrieve {
    label 'seqtk'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/illumina/", mode: 'copy', pattern: "*unmapped_*.fastq"
    input:
    tuple val(name), path(assembly), path(reads)
    output:
    path("unmapped_*.fastq") optionnal true
    shell:
    """
    ## illumina unmapped reads retrieval
    samtools view -f4 ${bam} > illumina_unmapped_contigs.sam
    cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
    seqtk subseq ${reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
    if [[ "${reads[1]}" != "" ]]; then
        seqtk subseq ${reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq
        fi

    rm *_unmapped_contigs.sam
    rm *_unmapped.list
    """

}


process unmapped_ont_retrieve {
    label 'seqtk'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/ont/", mode: 'copy', pattern: "*unmapped_*.fastq"
    input:
    tuple val(name), path(assembly), path(reads)
    output:
    path("unmapped_*.fastq") optionnal true
    shell:
    """
    ## ONT unmapped reads retrieval
    samtools view -f4 ${bam} > ont_unmapped_contigs.sam
    cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
    seqtk subseq ${reads} ont_unmapped.list > unmapped_ONT.fastq
    
    rm *_unmapped_contigs.sam
    rm *_unmapped.list
    """

}

// """
//     ## illumina unmapped reads retrieval
//     samtools view -f4 ${ill_bam} > illumina_unmapped_contigs.sam
//     cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
//     seqtk subseq ${ill_reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
//     seqtk subseq ${ill_reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq

//     ## ONT unmapped reads retrieval
//     samtools view -f4 ${ont_bam} > ont_unmapped_contigs.sam
//     cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
//     seqtk subseq ${ont_reads} ont_unmapped.list > unmapped_ONT.fastq

//     rm illumina_unmapped_contigs.sam
//     rm illumina_unmapped.list
//     rm ont_unmapped_contigs.sam
//     rm ont_unmapped.list
//     """