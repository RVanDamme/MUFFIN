process unicycler {
    label 'unicycler'
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.fasta"
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.gfa"
    errorStrategy { task.exitStatus in 1 ? 'finish' }
    input:
    set val(name),val(contig_list), file(illumina), file(ont)    
    output:
    set val(name), file("*_final.fasta") optionnal true
    file("*_final.gfa") optionnal true
    shell:
    """
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon
    mv output/assembly.fasta ${contig_list}"_final.fasta"
    mv output/assembly.gfa ${contig_list}"_final.gfa"
    """
}
