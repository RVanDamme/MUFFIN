process de_novo_transcript_and_quant {
    label 'trinity'
    // publishDir "${params.output}/${name}//", mode: 'copy', pattern: ""
    input:
    set val(name), file(rna)
    output:
    set val(name), file("trinity_transcript.fasta")
    set val(name), file("trinity_transcript_quant.sf")
    script:
    """
    mem=\$(echo !{task.memory} | sed 's/ GB/G/g')
    Trinity --seqType fq --max_memory \$mem --CPU !{task.cpus} --left !{rna[0]} --right !{rna[1]}
    align_and_estimate_abundance.pl --transcripts trinity_out_dir/Trinity.fasta --est_method salmon --left !{rna[0]} --right !{rna[1]} --seqType fq --output_dir quant_salmon --thread_count !{task.cpus}  --prep_reference
    cp trinity_out_dir/Trinity.fasta trinity_transcript.fasta
    cp quant_salmon/quant.sf trinity_transcript_quant.sf
    """
}

}