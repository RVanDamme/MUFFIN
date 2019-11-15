process unicycler {
    label 'unicycler'
    errorStrategy { task.exitStatus = 1 ? 'retry' : 'terminate' }
    maxRetries 2
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.fasta"
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.gfa"
    input:
    set val(name),val(contig_list), file(illumina), file(ont)    
    output:
    set val(name), file("*_final.fasta") optional true
    file("*_final.gfa") optional true
    shell:
    if (task.attempt == 1)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    if (task.attempt == 2)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp --max_kmer_frac 0.80
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    if (task.attempt == 3)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp --max_kmer_frac 0.60
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
}
