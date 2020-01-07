process racon {
    label 'racon'
    input:
        set val(name), file(read), file(assembly), file(mapping) 
    output:
        set val(name), file(read), file("${name}_consensus.fasta") 
    shell:
        """
        racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
        """
    }

process medaka {
    label 'medaka'
    input:
        set val(name), file(read), file(consensus) 
    output:
        set val(name), file("${name}_polished.fasta") 
    script:
        """
        medaka_consensus -i ${read} -d ${consensus} -o polished -t ${task.cpus} -m ${params.model}
        mv polished/consensus.fasta ${name}_polished.fasta
        """
  }

process pilon {
    label 'pilon'
    publishDir "${params.output}/${name}/flye_assembly/", mode: 'copy', pattern: "polished_assembly.fasta" 
    input:
        set val(name), file(assembly), file(ont_read)
        val(iteration)
    output:
        set val(name) , file("polished_assembly.fasta")
    shell:
    """
    mem=\$(echo !{task.memory} | sed 's/ GB//g')
    assemb="!{assembly}"
    for ite in {1..!{iteration}}
    do
        bwa index \$assemb
        bwa mem \$assemb !{ont_read} | samtools view -bS - | samtools sort -@ !{task.cpus} - > \$ite.bam
        samtools index -@ !{task.cpus} \$ite.bam
        pilon -Xmx\$mem"g" --threads !{task.cpus} --genome \$assemb --bam \$ite.bam --output \$ite"_polished_assembly"
        assemb=\$ite"_polished_assembly.fasta"
    done
    mv !{iteration}"_polished_assembly.fasta" polished_assembly.fasta
    """

}



//*********************************
// if polish with long on pilon add
//*********************************
// assemb2=\$bin_name"_illumina_polished.fasta"
//         for ite in {1..!{iteration}}
//         do
//             bwa index \$assemb2
//             bwa mem \$assemb2 !{ont_read} > assembly_ont_mapped.sam
//             samtools view -bS assembly_ont_mapped.sam > assembly_ont_mapped.bam
//             samtools sort -@ !{task.cpus} assembly_ont_mapped.bam > \$ite"_ont.bam"
//             samtools index -@ !{task.cpus} \$ite"_ont.bam"
//             pilon -Xmx24g --threads !{task.cpus} --genome \$assemb2 --bam \$ite"_ont.bam" --output \$ite"_polished_ont"
//             assemb=\$ite"_polished_ont.fasta"
//         done
//         mv !{iteration}"_polished_ont.fasta" \$bin_name"_ont_polished.fasta"


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