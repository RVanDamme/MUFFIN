process racon {
    label 'racon'

    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
        tuple val(name), path(read), path(assembly), path(mapping) 
    output:
        tuple val(name), path(read), path("${name}_consensus.fasta") 
    shell:
        """
        racon -t ${task.cpus} ${read} ${mapping} ${assembly} > ${name}_consensus.fasta
        """
    }
//medaka_consensus -i ${read} -d ${consensus} -o polished -t ${task.cpus} -m ${params.model}
process medaka {
    label 'medaka'

    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
        tuple val(name), path(read), path(consensus) 
    output:
        tuple val(name), path("${name}_polished.fasta") 
    script:
        """
        medaka --version
        medaka_consensus -i ${read} -d ${consensus} -o polished -t ${task.cpus} -m ${params.model}
        mv polished/consensus.fasta ${name}_polished.fasta
        """
  }

process pilon {
    label 'pilon'
    
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/pilon_polished/", mode: 'copy', pattern: "*polished_bin_assembly.fasta" 
    input:
        tuple val(name), path(assembly)
        tuple val(name), path(ill_read)
        val(iteration)
    output:
        tuple val(name) , path("*polished_bin_assembly.fasta")
    shell:
    """
    bin_id=\$(basename ${assembly} | sed -r "s/\\.\\w+//2")

    mem=\$(echo ${task.memory} | sed 's/ GB//g'| sed 's/g//g' | sed 's/ B//g')
    partial_mem=\$((\$mem*40/100))
    assemb="${assembly}"
    for ite in {1..${iteration}}
    do
        bwa index \$assemb
        bwa mem \$assemb ${ill_read[0]} ${ill_read[1]} | samtools view -bS - | samtools sort -@ ${task.cpus} - > \$ite.bam
        samtools index -@ ${task.cpus} \$ite.bam
        pilon -Xmx\$partial_mem"g" --threads ${task.cpus} --genome \$assemb --bam \$ite.bam --output \$ite"_polished_assembly"
        assemb=\$ite"_polished_assembly.fasta"
    done
    mv ${iteration}"_polished_assembly.fasta" \$bin_id"_polished_bin_assembly.fasta"
    """

}

process pilon2 {
    label 'pilon'
    
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/pilon_polished/", mode: 'copy', pattern: "pilon_res/*.fasta" 
    
    input:
        tuple val(name), path(assembly_files)
        tuple val(name), path(ill_read)
        val(iteration)
    
    output:
        tuple val(name), path("pilon_res/*.fasta")

    script:
    """
    mkdir -p ./pilon_res
    for assembly in ${assembly_files}
    do
        bin_id=\$(basename \$assembly | sed -r "s/\\.\\w+//2")
        mem=\$(echo ${task.memory} | sed 's/ GB//g'| sed 's/g//g' | sed 's/ B//g')
        partial_mem=\$((\$mem*40/100))
        assemb="\$assembly"
        for ite in {1..${iteration}}
        do
            bwa index \$assemb
            bwa mem \$assemb ${ill_read[0]} ${ill_read[1]} | samtools view -bS - | samtools sort -@ ${task.cpus} - > \$ite.bam
            samtools index -@ ${task.cpus} \$ite.bam
            pilon -Xmx\$partial_mem"g" --threads ${task.cpus} --genome \$assemb --bam \$ite.bam --output \$ite"_polished_assembly"
            assemb=\$ite"_polished_assembly.fasta"
        done
        mv ${iteration}"_polished_assembly.fasta" "./pilon_res/"\$bin_id"_polished_bin_assembly.fasta"
    done
    """
}

// process pilon_tk {
//     label 'pilon'

//     conda  'bioconda::pilon=1.23 bioconda::bwa=0.7.17 bioconda::samtools=1.9'
    
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     publishDir "${params.output}/${name}/assemble/assembly/pilon_polished/", mode: 'copy', pattern: "*polished_bin_assembly.fasta" 
//     input:
//         tuple val(name), path(assembly), path(ill_read)
//         val(iteration)
//     output:
//         tuple val(name) , path("polished_bin_assembly.fasta")
//     shell:
//     """
//     bin_id=\$(basename ${assembly} | sed -r "s/\\.\\w+//2")

//     bwa index -p illumina -a bwtsw ${assembly}
//     bwa mem illumina ${ill_read[0]} ${ill_read[1]} -t ${task.cpus} > illumina.sam
//     samtools view -bS illumina.sam > illumina.bam
//     samtools sort -@ ${task.cpus} -o illumina_sorted.bam illumina.bam

//     ## illumina mapped reads retrieval
//     samtools index -@ ${task.cpus} illumina_sorted.bam
//     samtools view -bh illumina_sorted.bam \$list > illumina_contigs.bam  
//     samtools view -F4 illumina_contigs.bam > illumina_mapped_contigs.sam
//     cut -f1 illumina_mapped_contigs.sam | sort | uniq > \$bin"_illumina_mapped.list"
//     seqtk subseq ${ill_reads[0]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R1.fastq"
//     seqtk subseq ${ill_reads[1]} \$bin"_illumina_mapped.list" > \$bin"_illumina_R2.fastq"

//     #do pilon
//     mem=\$(echo ${task.memory} | sed 's/ GB//g'| sed 's/g//g' | sed 's/ B//g')
//     partial_mem=\$((\$mem*40/100))
//     assemb="${assembly}"
//     for ite in {1..${iteration}}
//     do
//         bwa index \$assemb
//         bwa mem \$assemb \$bin"_illumina_R1.fastq" \$bin"_illumina_R2.fastq" | samtools view -bS - | samtools sort -@ ${task.cpus} - > \$ite.bam
//         samtools index -@ ${task.cpus} \$ite.bam
//         pilon -Xmx\$partial_mem"g" --threads ${task.cpus} --genome \$assemb --bam \$ite.bam --output \$ite"_polished_assembly"
//         assemb=\$ite"_polished_assembly.fasta"
//     done
//     mv ${iteration}"_polished_assembly.fasta" \$bin_id"_polished_bin_assembly.fasta"

//     #clean phase
//     rm illumina.*
//     rm illumina_contigs.bam
//     rm illumina_mapped_contigs.sam
//     rm *.bam.bai
//     """

// }

//use minimap instead of bwa
process pilong {
    label 'pilon'

    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/pilon_polished/", mode: 'copy', pattern: "polished_assembly.fasta" 
    input:
        tuple val(name), path(assembly), path(ont_read)
        val(iteration)
    output:
        tuple val(name) , path("polished_assembly.fasta")
    shell:
    """
    mem=\$(echo ${task.memory} | sed 's/ GB//g'| sed 's/g//g' | sed 's/ B//g')
    partial_mem=\$((\$mem*40/100))
    assemb="${assembly}"
    for ite in {1..${iteration}}
    do
        minimap2 -ax map-ont \$assemb ${ont_read} | samtools view -bS - | samtools sort -@ ${task.cpus} - > \$ite.bam
        samtools index -@ ${task.cpus} \$ite.bam
        pilon -Xmx\$partial_mem"g" --threads ${task.cpus} --genome \$assemb --bam \$ite.bam --output \$ite"_polished_assembly"
        assemb=\$ite"_polished_assembly.fasta"
    done
    mv ${iteration}"_polished_assembly.fasta" polished_assembly.fasta
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