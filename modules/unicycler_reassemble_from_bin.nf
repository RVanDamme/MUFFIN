process unicycler {
    maxForks 1
    label 'unicycler'
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*.fa"
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*.gfa"
    input:
    set val(name), val(bin_name), file(illumina), file(ont)    
    output:
    set val(name), val(bin_name), file("*.fa")
    file("*.gfa")
    shell:
    """
    mkdir spades_tmp
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon --spades_tmp_dir spades_tmp
    mv output/assembly.fasta ${bin_name}".fa"
    mv output/assembly.gfa ${bin_name}".gfa"
    """
}