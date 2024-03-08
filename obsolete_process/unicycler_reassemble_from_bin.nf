process unicycler {
    maxForks 1
    label 'unicycler'
    publishDir "${params.output}/${name}/assemble/reassembly/unicycler_bins/", mode: 'copy', pattern: "*.fa"
    publishDir "${params.output}/${name}/assemble/reassembly/unicycler_bins/", mode: 'copy', pattern: "*.gfa"
    errorStrategy { task.exitStatus in 1..1 ? 'retry' : 'finish'}
    maxRetries 4
    input:
    tuple val(name), val(bin_name), path(illumina), path(ont)    
    output:
    tuple val(name), path("*.fa") optional true
    path("*.gfa") optional true
    shell:
    if (task.attempt == 1)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    else if (task.attempt == 2)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp --max_kmer_frac 0.85
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    else if (task.attempt == 3)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp --max_kmer_frac 0.70
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    else if (task.attempt == 4)
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp --max_kmer_frac 0.50
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
    else
    error "Unicycler was unable to process your data please restart MUFFIN without Unicycler activated"
}