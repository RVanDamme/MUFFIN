#!/usr/bin/env nextflow

params.illumina = "./illumina"
params.reads_illumina = "${params.illumina}*_R{1,2}.fastq"
readFileList_illumina = Channel
                .fromFilePairs(params.reads_illumina)
                .ifEmpty { error "Cannot find any Illumina reads in the directory: ${params.reads} \n Delfault is ./illumina" }
params.nanopore ='./nanopore'
params.reads_nanopore = "${params.nanopore}*.fastq"
readFileList_nanopore = Channel
                .fromPath(params.reads_nanopore)
                .map {file -> tuple(file.baseName, file) }
                .ifEmpty { error "Cannot find any Nanopore reads in the directory: ${params.reads} \n Delfault is ./nanopore" }


params.polish = ''
params.output = ''
params.cpus = 4
params.mem = "8GB"
params.filter = FALSE
params.qual = FALSE

if (params.filter == '' & params.qual == '')
    readsFileList_ch = readFileList_illumina.join(readFileList_nanopore)
else if (params.filter == '' & params.qual == TRUE)
    readsFileList_qual = readFileList_illumina.join(readFileList_nanopore)
else if (params.filter == TRUE & params.qual == '')
    readsFileList_filter = readFileList_illumina.join(readFileList_nanopore)
else if (params.filter == TRUE & params.qual == TRUE)
    readsFileList_filter = readFileList_illumina.join(readFileList_nanopore)


// choose the Method (hybrid or long with short polishing)
params.assembler = ''
if (params.assembler != 'metaflye' && params.assembler != 'metaspades') {
        exit 1, "--method: ${params.method}. \
        Should be 'metaflye'or 'metaspades'"
}
// need --output
if (params.output == '') {
    exit 1, "--output is a required parameter"
}
if (params.polish == '') {
    println("No polishing; If you want polishing chose either 'medaka' 'racon' or 'both' in --polish")
}

//IF WE FILTER 
// if (params.filter == True) {}
//     process filter_reads {
//         input:
//         set val(name), file(illumina), file(nanopore) from readsFileList_filter

//         output:
//         if (params.qual == '')
//              set val(name), file("${name}_filtered_1.fasta"), file("${name}_filtered_2.fasta"), file("${name}_filtered_nanopore.fasta") into {readsFileList_qual}
//         if (params.qual == TRUE)
//              set val(name), file("${name}_filtered_1.fasta"), file("${name}_filtered_2.fasta"), file("${name}_filtered_nanopore.fasta") into {readsFileList_ch}
//         script:
//             """
//             NEED TO LOOK FOR IT
//              map using bbmap or bowtie and only keep the unmap
//              use blast and only keep the non hit
//             """
//     }

//IF WE Qual FILTER 
// if (params.qual == TRUE) {}
//     process filter_reads {
//         input:
//         set val(name), file(illumina), file(nanopore) from readsFileList_qual

//         output:
//         set val(name), file("${name}_filtered_1.fasta"), file("${name}_filtered_2.fasta"), file("${name}_filtered_nanopore.fasta") into {readsFileList_ch}
//   BIG ISSUES THE ILLUMINA ARE  A TUPLE AND USE AS A TUPLE EVERYWHERE SO NEED TO PUT THEM AS A TUPLE HERE ALSO
//    MAYBE USING file("${name}_filtered_1.fasta","${name}_filtered_2.fasta") will WORK need to test it
//         script:
//             """
//             NEED TO LOOK FOR IT
//                for nanopore use FiltLong
//                for illumina use fastp?

//             """
//     }



process assembly {
    echo true
    input:
        set val(name), file(illumina), file(nanopore) from readsFileList_ch  
    output:
        set val(name), file(illumina), file(nanopore), file('assembly.fasta') into {qc_assembly, mapping_assembly, polishing_assembly }

    script:
    if(params.assembler == 'metaspades')
        """
        spades.py -1 ${illumina[0]} -2 ${illumina[1]}  --meta --nanopore ${nanopore} -o spades_output -t ${task.cpus} -m ${task.mem}
        mv spades_output/contigs.fasta  assembly.fasta

        """
    else if(params.assembler == 'metaflye')
        """
        flye --nano-corr ${nanopore} -o flye_output -t ${task.cpus} --plasmids --meta
        mv flye_output/assembly.fasta assembly.fasta

        """
        //NEED TO ADD THE POLISHING STEP (one round Medaka then Pilon or racon)

}  

if (params.assembler == 'metaflye')
    process polishing {
        input:
            set val(name), file(illumina), file(nanopore), file(assembly) from polishing_assembly

        output:
            set val(name), file(illumina), file(nanopore), file('assembly_polished.fasta') into {qc_polished, mapping_polished}

        script:
        if (params.polish == 'medaka')
            """
            medaka_consensus -i ${nanopore} -d ${assembly} -o consensus -t ${task.cpus}
            mv consensus/consensus.fasta medaka_consensus.fasta
            """
        else if (params.polish == 'racon')
            """
            minimap2 -ax map-ont ${assembly} ${nanopore} > ready_polish.sam
            racon ${nanopore} ready_polish.sam ${assembly}

            """
        else if (params.polish == 'both')
            """
            minimap2 -ax map-ont ${assembly} ${nanopore} > ready_to_polish.sam
            racon ${nanopore} ready_polish.sam ${assembly} -t ${task.cpus} > assembly_racon.fasta
            medaka_consensus -i ${nanopore} -d assembly_racon.fasta -o consensus -t ${task.cpus}
            mv consensus/consensus.fasta polished_consensus.fasta

            """
    }

// process assembly_qc {
//     input:
//         if (params.polish == "")
//             set val(name), file(illumina), file(nanopore), file(assembly) from qc_assembly      
//         else 
//             set val(name), file(illumina), file(nanopore), file(assembly) from qc_polished
//     script:
//     """
//     python metaquast.py ${assembly} 

//     """
// }

process assembly_mapping {
    echo true
    input:
        if (params.polish == "")
            set val(name), file(illumina), file(nanopore), files(assembly) from mapping_assembly
        else
            set val(name), file(illumina), file(nanopore), files(assembly) from mapping_polished
    output:
        set val(name), file(illumina), file(nanopore), files(assembly), files("mapping_illumina.sam"), files("mapping_nanopore.sam") into binning_assembly

    script:
        """
        minimap2 -ax map-ont ${assembly} ${nanopore} > mapping_nanopore.sam
        
        bwa index -p illumina -a bwtsw ${assembly}
        bwa mem illumina ${illumina} -t ${task.cpus} > mapping_illumina.sam

        """
}



process binning {
    echo true
    input:
    set val(name), file(illumina), file(nanopore), files(assembly), files(mapping_illumina), files(mapping_nanopore) from binning_assembly   
        
    output:
        

    script:
    """
    metabat2  -o output_bin -t ${task.cpus} --unbinned --seed 133742 -i ${assembly} ${mapping_illumina} ${mapping_nanopore}

    """ // NEED TO FIND HOW TO ADD MAPPING AND IF RELEVANT
}

