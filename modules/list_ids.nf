process contig_list {
    label ''
    input:
    set val(name), val(bins)
    output:
    set val(name), file("ids/*.list")
    script:
    """
    mkdir ids
    for bin in \$(ls -d !{refined_bins}/bin*);
        do
        cat \$bin/*.{fasta;fa} | grep -o -E "^>\w+" | tr -d "@" > "ids/"\$bin".contigs.list" ;
        done ;
    """
}

process reads_retrieval {
    label ''
    input:
    set val(name),file(contig_list), file(ill_bam), file(ont_bam), file(ill_reads), file(ont_reads)
    output:
    set val(name), val("\$bin"), file(["\$bin_illumina_R{1,2}.fastq"]), file("\$bin_ont.fastq")
    script:
    // first I extract the reads that NEED TO REDO IT WITH FRESH MIND (include BWA.nf)
    """
    bin=\$(basename !{contig_list})
    list=\$(cat !{contig_list} | tr "\n" " " )
    
    ## illumina reads retrieval
    samtools view -bh !{ill_bam} \$list > ont_contigs.bam  
    samtools view -F4 ont_contigs.bam > ont_mapped_contigs.sam
    cut -f1 ill_mapped_contigs.sam | sort | uniq > \$bin"_illumina_mapped.list"
    seqtk subseq !{illumina[0]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R1.fastq"
    seqtk subseq !{illumina[1]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R2.fastq"

    ## ONT reads retrieval
    samtools view -bh !{ont_bam} \$list > ont_contigs.bam  
    samtools view -F4 ont_contigs.bam > ont_mapped_contigs.sam
    cut -f1 ont_mapped_contigs_sam | sort | uniq > \$bin"_ont_mapped.list"
    seqtk subseq !{illumina[0]} \$bin"_illumina_mapped.list" > \$bin"_ont.fastq"

    """
}

// notes that contigs.bam is the bam file of the reads aligned to the list of contigs used
// notes that mapped.contigs.bam is the mapped reads to the contigs