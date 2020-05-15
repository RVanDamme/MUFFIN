process spades {
    label 'spades'
    publishDir "${params.output}/${name}/spades_assembly/", mode: 'copy', pattern: "assembly.fasta" 
    input:
    tuple val(name), file(illumina), file(ont)
    output:
    tuple val(name), file("assembly.fasta")
    
    script:
    """
    mem=\$(echo !{task.memory} | sed 's/g//g')
    spades.py -1 !{illumina[0]} -2 !{illumina[1]}  --meta --nanopore !{ont} -o spades_output -t !{task.cpus} -m \$mem
    mv spades_output/contigs.fasta  assembly.fasta
    """

}