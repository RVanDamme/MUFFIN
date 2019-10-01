process spades {
    label 'spades'
    publishDir "${params.output}/${name}_assembly/", mode: 'copy', pattern: "assembly.fasta"
    input:
    set val(name), file(illumina), file(ont)
    output:
    set val(name), file("assembly.fasta")
    
    script:
    """
    spades.py -1 ${illumina[0]} -2 ${illumina[1]}  --meta --nanopore ${ont} -o spades_output -t ${task.cpus} -m ${task.memory}
    mv spades_output/contigs.fasta  assembly.fasta
    """
}