process flye {
    label 'flye'
    label 'ubuntu'
    if (params.assembly == true ) {publishDir "${params.output}/${name}_assembly/", mode: 'copy', pattern: "assembly.fasta"}
    input:
    set val(name), file(ont), file(genome_size)
    output:
    set val(name), file("assembly.fasta")
    shell:
    """
    size=\$(cat !{genome_size})
    flye --nano-corr !{nanopore} -o flye_output -t !{task.cpus} --plasmids --meta --genome-size \$size
    mv flye_output/assembly.fasta assembly.fasta
    """

}