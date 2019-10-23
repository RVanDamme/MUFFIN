process reads_retrieval {
    label 'ubuntu'
    label 'seqtk'
    if (params.out_bin_reads == true ) {publishDir "${params.output}/${name}_bins_reads/", mode: 'copy', pattern: "\$bin*.fastq"}
    if (params.out_unmapped == true ) {publishDir "${params.output}/${name}_unmapped_bam/", mode: 'copy', pattern: "*unmapped_*.fastq"}
    input:
    set val(name), file(contig_list), file(ill_bam), file(ont_bam), file(ill_reads), file(ont_reads), file(assembly)
    output:
    set val(name), file(contig_list), file("*_illumina_R{1,2}.fastq"), file("*_ont.fastq")
    shell:
    // first I extract the reads that NEED TO REDO IT WITH FRESH MIND (include BWA.nf)
    """
    bin=\$(basename !{contig_list})
    list=\$(cat !{contig_list} | tr "\n" " " ) 
    
    ## illumina mapped reads retrieval
    samtools index !{ill_bam}
    samtools view -bh !{ill_bam} \$list > illumina_contigs.bam  
    samtools view -F4 illumina_contigs.bam > illumina_mapped_contigs.sam
    cut -f1 illumina_mapped_contigs.sam | sort | uniq > \$bin"_illumina_mapped.list"
    seqtk subseq !{ill_reads[0]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R1.fastq"
    seqtk subseq !{ill_reads[1]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R2.fastq"

    ## illumina unmapped reads retrieval
    samtools view -f4 !{ill_bam} > illumina_unmapped_contigs.sam
    cut -f1 illumina_unmapped_contigs.sam | sort | uniq > illumina_unmapped.list
    seqtk subseq !{ill_reads[0]} illumina_unmapped.list > unmapped_ILL_R1.fastq
    seqtk subseq !{ill_reads[1]} illumina_unmapped.list > unmapped_ILL_R2.fastq

    ## ONT mapped reads retrieval
    samtools index !{ont_bam}
    samtools view -bh !{ont_bam} \$list > ont_contigs.bam  
    samtools view -F4 ont_contigs.bam > ont_mapped_contigs.sam
    cut -f1 ont_mapped_contigs.sam | sort | uniq > \$bin"_ont_mapped.list"
    seqtk subseq !{ont_reads} \$bin"_ont_mapped.list" > \$bin"_ont.fastq"

    ## ONT unmapped reads retrieval
    samtools view -f4 !{ont_bam} > ont_unmapped_contigs.sam
    cut -f1 ont_unmapped_contigs.sam | sort | uniq > ont_unmapped.list
    seqtk subseq !{ont_reads} ont_unmapped.list > unmapped_ONT.fastq
    """
}

// notes that contigs.bam is the bam file of the reads aligned to the list of contigs used
// notes that mapped.contigs.bam is the mapped reads to the contigs