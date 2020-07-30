process spades {
    label 'spades'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/spades/", mode: 'copy', pattern: "assembly.fasta" 
    input:
    tuple val(name), path(illumina), path(ont)
    output:
    tuple val(name), path("assembly.fasta")
    
    script:
    """
    mem=\$(echo ${task.memory} |sed 's/ GB//g' | sed 's/g//g')
    cpus=\$(echo ${task.cpus})
    echo \$cpus \$mem
    spades.py -1 ${illumina[0]} -2 ${illumina[1]}  --meta --nanopore ${ont} -o spades_output -t \$cpus -m \$mem
    mv spades_output/contigs.fasta  assembly.fasta
    """

}