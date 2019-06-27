#!/usr/bin/env nextflow

params.reads_illumina_1 = ''
params.reads_illumina_2 = ''
params.reads_nanopore = ''

params.output = ''
params.cpus = 4
params.mem = "8GB"

// choose the Method (hybrid or long with short polishing)
params.method = ''
if (params.method != 'meta-flye' && params.assembler != 'meta-spades') {
        exit 1, "--method: ${params.method}. \
        Should be 'meta-flye'or 'meta-spades'"
}

// need --output
if (params.output == '') {
    exit 1, "--output is a required parameter"
}
//requires --illumina_1
if (params.reads_illumina_1 == '') {
    exit 1, "--illumina_1 is a required parameter"
}
//requires --illumina_2
if (params.reads_illumina_2 == '') {
    exit 1, "--illumina_2 is a required parameter"
}
//requires --nanopore
if (params.reads_nanopore == '') {
    exit 1, "--nanopore is a required parameter"
}

// process filter_illumina {
//     input:
// 	   file(reads_illumina_1) from file(params.reads_illumina_1)
// 	   file(reads_illumina_2) from file(params.reads_illumina_2)
//     output:
// 	   file('filtered_illumina_1.fastq') into filtered_reads_illumina_1
// 	   file('filtered_illumina_2.fastq') into filtered_reads_illumina_2

// 	script:
//         """
//     	NEED TO LOOK FOR IT
//         """
// }


// process filter_nanopore {
//     input:
//     file(reads_nanopore) from file(params.reads_nanopore)
//     output:
// 	   file('filtered_nanopore.fastq') into filtered_reads_nanopore
// 	script:
//         """
//     	NEED TO LOOK FOR IT
//         """
// }


// need the reads in assembly, and mapping so the channel is dupe
// filtered_reads_illumina_1.into { filtered_reads_illumina_assembly_1; filtered_reads_illumina_mapping_1 }
// filtered_reads_illumina_2.into { filtered_reads_illumina_assembly_2; filtered_reads_illumina_mapping_2 }
// filtered_reads_nanopore.into { filtered_reads_nanopore_assembly; filtered_reads_nanopore_mapping }


// need the reads in assembly, and mapping so the channel is dupe
// as long as we don't filter we use the "raw" reads
reads_illumina_1.into { filtered_reads_illumina_assembly_1; filtered_reads_illumina_mapping_1 }
reads_illumina_2.into { filtered_reads_illumina_assembly_2; filtered_reads_illumina_mapping_2 }
reads_nanopore.into { filtered_reads_nanopore_assembly; filtered_reads_nanopore_mapping }

process assembly {
    //publishDir params.output, mode: 'copy', pattern: 'assembly.fasta'  ?????

    input:
        file reads_illumina_1 from filtered_reads_illumina_assembly_1
        file reads_illumina_2 from filtered_reads_illumina_assembly_2
        file reads_nanopore from filtered_reads_nanopore_assembly
        
    output:
        file 'assembly.fasta' into assembly

    script:
    if(params.method == 'meta-spades')
        """
        spades.py -1 ${reads_illumina_1} -2 ${reads_illumina_2}  --meta --nanopore ${reads_nanopore} -o spades_output -t ${task.cpus} -m ${task.mem}
        mv spades_output/contigs.fasta  assembly.fasta

        """
    else if(params.assembler == 'meta-flye')
        """
        flye --nano-corr ${reads_nanopore} -o flye_output -t ${task.cpus} --plasmids --meta
        mv flye_output/assembly.fasta assembly.fasta

        """
        //NEED TO ADD THE POLISHING STEP (probably racon)

}

assembly.into {qc_assembly; mapping_assembly; binning_assembly}

// process assembly_qc {
//     input:
//         file assembly from qc_assembly
        
        
//     output:
        

//     script:
//     """
//     TO BE DETERMINED

//     """

// }

process assembly_mapping {
    input:
        file assembly from mapping_assembly
        file reads_illumina_1 from filtered_reads_illumina_mapping_1
        file reads_illumina_2 from filtered_reads_illumina_mapping_2
        file reads_nanopore from filtered_reads_nanopore_mapping
        
    output:
        file 'mapping_illumina.bam' into mapping_illumina
        file 'mapping_nanopore.bam' into mapping_nanopore

    script:
    """


    """
}

process binning {
    input:
        file assembly from binning_assembly
        file mapping_illumina from mapping_illumina
        file mapping_nanopore from mapping_nanopore
        
    output:
        

    script:
    """
    metabat2 -i ${assembly} -o output_bin -t ${task.cpus} --unbinned --seed 133742

    """ // NEED TO FIND HOW TO ADD MAPPING AND IF RELEVANT
}