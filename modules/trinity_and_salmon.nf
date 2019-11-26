process de_novo_transcript {
    label 'trinity'
    // publishDir "${params.output}/${name}//", mode: 'copy', pattern: ""
    input:
    set val(name), file(rna)
    output:
    set val(name), file("trinity_transcript.fasta")
    shell:
    """
    mem=\$(echo "!{task.memory}" | sed 's/ GB/G/g')
    echo \$mem
    Trinity --seqType fq --max_memory 20G --CPU !{task.cpus} --left !{rna[0]} --right !{rna[1]}
    cp trinity_out_dir/Trinity.fasta trinity_transcript.fasta
    """
}
process quantification {
    label 'trinity'
    // publishDir "${params.output}/${name}//", mode: 'copy', pattern: ""
    input:
    set val(name), file(rna)
    set val(name), file(dammit_fa), file(dammit_gff)
    output:
    set val(name), file(dammit_fa), file(dammit_gff), file("trinity_transcript_quant.sf")
    shell:
    """
    mem=\$(echo "!{task.memory}" | sed 's/ GB/G/g')
    echo \$mem
    align_and_estimate_abundance.pl --transcripts !{dammit_fa} --est_method salmon --left !{rna[0]} --right !{rna[1]} --seqType fq --output_dir quant_salmon --thread_count !{task.cpus}  --prep_reference
    cp quant_salmon/quant.sf trinity_transcript_quant.sf
    """
}
