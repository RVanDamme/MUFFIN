process unicycler {
    label 'unicycler'
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.fasta"
    publishDir "${params.output}/${name}/unicycler_assembly/", mode: 'copy', pattern: "*_final.gfa"
    errorStrategy { task.exitStatus = 1 ? 'finish' : 'terminate' }
    input:
    set val(name),val(contig_list), file(illumina), file(ont)    
    output:
    set val(name), file("*_final.fasta") optional true
    file("*_final.gfa") optional true
    shell:
    """
    unicycler -1 ${illumina[0]} -2 ${illumina[1]} -l ${ont} -o output -t ${task.cpus} --keep 0 --no_pilon
    mv output/assembly.fasta ${contig_list}"_final.fasta"
    mv output/assembly.gfa ${contig_list}"_final.gfa"
    """
}
