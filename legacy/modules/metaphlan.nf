

process metaphlan {
    label 'metaphlan' 

    publishDir "${params.output}/${name}/classify/metaphlan/", mode: 'copy', pattern: "*.txt"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(bin)
    output:
    path('*.txt')
    shell:
    """
    bin_id=\$(basename ${bin} | sed -r "s/\\.\\w+//2")
    ${bins} --input_type fasta --nproc ${task.cpus} > \$bin_id.txt
    """
}

