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


readFileList_ch = readFileList_illumina.join(readFileList_nanopore)

params.polish = ''
params.output = ''
params.cpus = 4
params.mem = "8GB"

// choose the Method (hybrid or long with short polishing)
params.method = ''
if (params.method != 'metaflye' && params.assembler != 'metaspades') {
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
// if (params.filter == '') {}
//     process filter_reads {
//         input:
//         set val(name), file(illumina), file(nanopore) from readsFileList_ch

//         output:
//         set val(name), file("${name}_filtered_1.fasta"), file("${name}_filtered_2.fasta"), file("${name}_filtered_nanopore.fasta") into {reads_assembly}

//         script:
//             """
//             NEED TO LOOK FOR IT
//             """
//     }




process assembly {
    echo true
    input:
        set val(name), file(illumina), file(nanopore) from readsFileList_ch  
    output:
        set val(name), file(illumina), file(nanopore), file('assembly.fasta') into {qc_assembly, mapping_assembly, polishing_assembly }

    script:
    if(params.method == 'metaspades')
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
            medaka_consensus -i ${nanopore} -d ${assembly} -o consensus -t ${task.cpus}
            mv consensus/consensus.fasta medaka_consensus.fasta
            minimap2 -ax map-ont medaka_consensus.fasta ${nanopore} > ready_to_polish.sam
            racon ${nanopore} ready_polish.sam medaka_consensus.fasta -t ${task.cpus} > assembly_polished.fasta
            """
    }

process assembly_qc {
    input:
        if (params.polish == "")
            set val(name), file(illumina), file(nanopore), file(assembly) from qc_assembly      
        else 
            set val(name), file(illumina), file(nanopore), file(assembly) from qc_polished
    script:
    """
    python metaquast.py ${assembly} 

    """
}

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

