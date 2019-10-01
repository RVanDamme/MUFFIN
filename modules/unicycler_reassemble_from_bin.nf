process unicycler {
    label 'unicycler'
    publishDir "${params.output}/${name}_contig_reassemblies_from_bins/", mode: 'copy', pattern: "${bin_name}*.{fasta,gfa}"
    input:
    set val(name), val(bin_name), file(illumina), file(ont)    
    output:
    set val(name), file("${bin_name}_polished.fasta"), file("${bin_name}.gfa")
    script:
    """
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o ${bin_name}_output -t ${task.cpus} --keep 0
    mv ${bin_name}_output/assembly.fasta ${bin_name}_polished.fasta
    mv ${bin_name}_output/assembly.gfa ${bin_name}.gfa
    """
}