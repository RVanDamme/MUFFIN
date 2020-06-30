process reads_retrieval {
    label 'seqtk'
    publishDir "${params.output}/${name}/assembled/reassembly/mapped_reads/", mode: 'copy', pattern: "*.fastq"
    input:
    tuple val(name), path(contig_list), path(ill_bam), path(ont_bam), path(ill_reads), path(ont_reads)
    output:
    tuple val(name), val(file(file(file(contig_list).baseName).baseName).baseName), path("*_illumina_R{1,2}.fastq"), path("*_ont.fastq")
    shell:
    // first I extract the reads that NEED TO REDO IT WITH FRESH MIND (include BWA.nf)
    """
    bin=\$(basename -s .fa.contigs.list ${contig_list})
    list=\$(cat ${contig_list} | tr "\\n" " " ) 

    
    ## illumina mapped reads retrieval
    samtools index -@ ${task.cpus} ${ill_bam}
    samtools view -bh ${ill_bam} \$list > illumina_contigs.bam  
    samtools view -F4 illumina_contigs.bam > illumina_mapped_contigs.sam
    cut -f1 illumina_mapped_contigs.sam | sort | uniq > \$bin"_illumina_mapped.list"
    seqtk subseq ${ill_reads[0]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R1.fastq"
    seqtk subseq ${ill_reads[1]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R2.fastq"

    ## ONT mapped reads retrieval
    samtools index -@ ${task.cpus} ${ont_bam}
    samtools view -bh ${ont_bam} \$list > ont_contigs.bam  
    samtools view -F4 ont_contigs.bam > ont_mapped_contigs.sam
    cut -f1 ont_mapped_contigs.sam | sort | uniq > \$bin"_ont_mapped.list"
    seqtk subseq ${ont_reads} \$bin"_ont_mapped.list" > \$bin"_ont.fastq"

    rm illumina_contigs.bam
    rm illumina_mapped_contigs.sam
    rm ont_contigs.bam
    rm ont_mapped_contigs.sam
    rm *.bam.bai
    """

}

// TODO PUT THE UNMAPPING AS A STAND ALONE STEP RETRIEVEING ALL UNMAPPED AS ONE FILE

// notes that contigs.bam is the bam file of the reads aligned to the list of contigs used
// notes that mapped.contigs.bam is the mapped reads to the contigs

process unmapped_retrieve {
    label 'seqtk'
    publishDir "${params.output}/${name}/assembled/reassembly/unmapped_reads/", mode: 'copy', pattern: "*unmapped_*.fastq"
    input:
    tuple val(name), path(ill_bam), path(ont_bam), path(ill_reads), path(ont_reads)
    output:
    path("unmapped_*.fastq") optionnal true
    shell:
    """
    ## illumina unmapped reads retrieval
    samtools view -f4 ${ill_bam} > illumina_unmapped_contigs.sam
    cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
    seqtk subseq ${ill_reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
    seqtk subseq ${ill_reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq

    ## ONT unmapped reads retrieval
    samtools view -f4 ${ont_bam} > ont_unmapped_contigs.sam
    cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
    seqtk subseq ${ont_reads} ont_unmapped.list > unmapped_ONT.fastq

    rm illumina_unmapped_contigs.sam
    rm illumina_unmapped.list
    rm ont_unmapped_contigs.sam
    rm ont_unmapped.list
    """

}