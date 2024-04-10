// process unmapped_retrieve {
//     label 'seqtk'
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/", mode: 'copy', pattern: "*unmapped_*.fastq"
//     input:
//     tuple val(name), path(bam), path(reads), val(sequencing_type)
//     output:
//     path("unmapped_*.fastq") optionnal true
//     shell:
//     """
//     if [[ "${sequencing_type}" == "illumina" ]]; then
//         ## illumina unmapped reads retrieval
//         samtools view -f4 ${bam} > illumina_unmapped_contigs.sam
//         cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
//         seqtk subseq ${reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
//         if [[ "${reads[1]}" != "" ]]; then
//             seqtk subseq ${reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq
//             fi
       
//     elif [[ "${sequencing_type}" == "ont" ]]; then
//         ## ONT unmapped reads retrieval
//         samtools view -f4 ${bam} > ont_unmapped_contigs.sam
//         cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
//         seqtk subseq ${reads} ont_unmapped.list > unmapped_ONT.fastq
//     fi

//     rm *_unmapped_contigs.sam
//     rm *_unmapped.list
//     """

// }

process unmapped_illumina_retrieve {
    label 'seqtk'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/illumina/", mode: 'copy', pattern: "*unmapped_*.fastq"
    input:
    tuple val(name), path(bam), path(reads)
    output:
    path("unmapped_*.fastq")
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
    tuple val(name), path(bam), path(reads)
    output:
    path("unmapped_*.fastq")
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

process ont_reads_retrieval {
    label 'seqtk'
    publishDir "${params.output}/${name}/assembled/reassembly/mapped_reads/ont/", mode: 'copy', pattern: "*.fastq"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(ont_bam), path(ont_reads)
    output:
    tuple val(name), path("*_ont.fastq")
    shell:
    // first I extract the reads that NEED TO REDO IT WITH FRESH MIND (include BWA.nf)
    """
    ## ONT mapped reads retrieval
    samtools index -@ ${task.cpus} ${ont_bam}
    samtools view -bh ${ont_bam} \$list > ont_contigs.bam  
    samtools view -F4 ont_contigs.bam > ont_mapped_contigs.sam
    cut -f1 ont_mapped_contigs.sam | sort | uniq > \$bin"_ont_mapped.list"
    seqtk subseq ${ont_reads} \$bin"_ont_mapped.list" > \$bin"_ont.fastq"

    rm ont_contigs.bam
    rm ont_mapped_contigs.sam
    rm *.bam.bai
    """

}

process illumina_reads_retrieval {
    label 'seqtk'
    publishDir "${params.output}/${name}/assembled/reassembly/mapped_reads/illumina/", mode: 'copy', pattern: "*.fastq"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(illumina_bam), path(illumina_reads)
    output:
    tuple val(name), path("*_ont.fastq")
    shell:
    // first I extract the reads that NEED TO REDO IT WITH FRESH MIND (include BWA.nf)
    """
    ## illumina mapped reads retrieval
    samtools index -@ ${task.cpus} ${ill_bam}
    samtools view -bh ${ill_bam} \$list > illumina_contigs.bam  
    samtools view -F4 illumina_contigs.bam > illumina_mapped_contigs.sam
    cut -f1 illumina_mapped_contigs.sam | sort | uniq > \$bin"_illumina_mapped.list"
    seqtk subseq ${ill_reads[0]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R1.fastq"
    seqtk subseq ${ill_reads[1]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R2.fastq"

    rm illumina_contigs.bam
    rm illumina_mapped_contigs.sam
    rm *.bam.bai
    """

}