process de_novo_transcript_and_quant {
    maxForks 1
    label 'trinity'

    publishDir "${params.output}/${name}/annotate/de_novo_transcript/", mode: 'copy', pattern: "*_transcript.fasta"
    publishDir "${params.output}/${name}/annotate/quant_of_transcript/", mode: 'copy', pattern: "*_transcript_quant.sf"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(rna)
    output:
    tuple val(name), path("*_transcript.fasta"), path("*_transcript_quant.sf")
    shell:
    """
    mem=\$(echo "!{task.memory}" | sed 's/ GB/g/g' | sed 's/g/G/g' | sed 's/ B/G/g')
    echo \$mem
    Trinity --seqType fq --max_memory \$mem --CPU !{task.cpus} --left !{rna[0]} --right !{rna[1]}
    cp trinity_out_dir/Trinity.fasta !{name}_transcript.fasta
    align_and_estimate_abundance.pl --transcripts !{name}_transcript.fasta --est_method salmon --left !{rna[0]} --right !{rna[1]} --seqType fq --output_dir quant_salmon --thread_count !{task.cpus}  --prep_reference
    cp quant_salmon/quant.sf !{name}_transcript_quant.sf
    """
}
