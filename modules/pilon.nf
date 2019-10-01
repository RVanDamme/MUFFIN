process mapper {
    label 'pilon'
    label 'ubuntu'
    publishDir "${params.output}/tmp/${name}_assembly/", mode: 'copy', pattern: "polished_assembly.fasta"
    input:
    set val(name), file(ont_bam)
    set val(name), file(assembly)
    set val(iteration)
    set val(threshold) // which track (--tracks) to look at as a threshold?
    output:
    set val(name) , file(polished_assembly.fasta)
    shell:
    """
    assemb="!{assembly}"
    for ite in {1..!{iteration}}
    do
        pilon --genome \$assemb --bam !{ont_bam} --output \$ite"_polished_assembly"
        assemb=\$ite"_polished_assembly.fasta"
    done
    mv \$iteration"_polished_assembly.fasta" polished_assembly.fasta
    """
}
    // *************
    // for threshold
    // *************

    // assemb="!{assembly}"
    // if !{threshold}==""
    // then
    //     for ite in {1..!{iteration}}
    //     do
    //          pilon --genome \$assemb --bam !{ont_bam} --output \$ite"_polished_assembly"
    //          assemb=\$ite"_polished_assembly.fasta"
    //     done
    //     mv \$iteration"_polished_assembly.fasta" polished_assembly.fasta
    // fi

    // if !{threshold}!=""
    // then
    //     thresh=\$(VALUE)
    //     while \$thresh -lt !{threshold}
    //          do
    //          pilon --genome \$assemb --bam !{ont_bam} --output \$thresh"_polished_assembly"
    //          assemb=\$thresh"_polished_assembly.fasta"
    //          thresh=\$(GET NEW VAL FROM \$thresh"_polished_assembly")
    //          done
    // fi