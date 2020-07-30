process flye {
    label 'flye'
    publishDir "${params.output}/${name}/assemble/assembly/flye_unpolished", mode: 'copy', pattern: "assembly.fasta"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(ont), path(genome_size)
    output:
    tuple val(name), path("assembly.fasta")
    shell:
    """
    size=\$(cat !{genome_size})
    flye --nano-raw ${ont} -o flye_output -t ${task.cpus} --plasmids --meta --genome-size \$size
    mv flye_output/assembly.fasta assembly.fasta
    """

}

//for flye updated over 2.7 use --nano-raw