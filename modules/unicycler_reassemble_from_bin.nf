process unicycler {
    label 'unicycler'
    publishDir "${params.output}/${name}_contig_unicycler_from_bins/", mode: 'copy', pattern: "${bin_name}*.{fasta,gfa}"
    input:
    set val(name),file(contig_list), file(illumina), file(ont)    
    output:
    set val(name), file("*_final.fasta"), file(illumina), file(ont)
    shell:
    """
    bin_name=\$(basename !{ont} | sed 's/_ont.fastq//g')
    unicycler -1 !{illumina[0]} -2 !{illumina[1]} -l !{ont} -o output -t !{task.cpus} --keep 0 --no_pilon
    mv output/assembly.fasta \$bin_name"_final.fasta"
    mv output/assembly.gfa \$bin_name"_final.gfa"
    """
}