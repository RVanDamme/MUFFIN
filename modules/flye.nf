process flye {
    label 'flye'
    publishDir "${params.output}/${name}/flye_assembly/", mode: 'copy', pattern: "assembly.fasta"
    input:
    tuple val(name), file(ont), file(genome_size)
    output:
    tuple val(name), file("assembly.fasta")
    shell:
    """
    size=\$(cat !{genome_size})
    flye --nano-raw ${ont} -o flye_output -t ${task.cpus} --plasmids --meta --genome-size \$size
    mv flye_output/assembly.fasta assembly.fasta
    """

}

//for flye updated over 2.7 use --nano-raw