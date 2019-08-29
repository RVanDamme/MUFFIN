#!/usr/bin/env nextflow

// TODO IMPLEMENT OUTPUT; FILTER ADN QC IF WANTED
//-with-dag file.html

params.assembler = ''
params.polish = "skip"
params.output = './'
//params.cpus = 4
//params.mem = "16GB"
params.genome = ""

// params.filter = FALSE
// params.qual = FALSE

if (params.assembler!='metaflye' && params.assembler!='metaspades') {
    exit 1, "--assembler: ${params.method}. Should be 'metaflye' or 'metaspades'"
}

if (params.assembler=='metaspades' && params.polish!="skip") {
    exit 1, "Do not specify '--polish' if you use '--assembler metaspades'"
}

if (params.assembler=='metaflye' && params.genome==""){
    exit 1, "please specify an estimate genome size using --genome and value (e.g. 5m or 2.6g) \n the genome size is required but metaspades is not very sensitive to it"
}

if (params.assembler=='metaflye' && params.polish=="skip") {
    println "No polishing step, if you want one use --polish"
}

if (params.assembler=='metaflye' && params.polish==true) {
    exit 1, "If you want to polish please specify the method with:\n --polish racon\n --polish medaka\n --polish both"
}

if (params.polish!="skip" && params.polish!='medaka' && params.polish!='racon' && params.polish!='both'){
    exit 1, "If you want to polish please specify the method with one of this:\n --polish racon\n --polish medaka\n --polish both"
}


params.filter = false
params.qual = false
params.illumina = "./data/tmft"
params.reads_illumina = "${params.illumina}*_R{1,2}.{fastq,fastq.gz}"
readFileList_illumina = Channel.fromFilePairs(params.reads_illumina).ifEmpty { error "Cannot find any Illumina reads in the directory: ${params.illumina} \n Delfault is ./illumina" }
params.nanopore ='./nanopore'
//params.reads_nanopore = "${params.nanopore}*.fastq"
params.reads_nanopore= "${params.nanopore}*.{fastq,fastq.gz}"
readFileList_nanopore = Channel.fromPath(params.reads_nanopore).map {file -> tuple(file.baseName, file) }.ifEmpty { error "Cannot find any Nanopore reads in the directory: ${params.nanopore} \n Delfault is ./nanopore" }

readFileList_nanopore.into {readFileList_nanopore; test}
test.subscribe{println "got: ${it}"}


if (params.filter == false & params.qual == false){
    readsFileList2_ch = readFileList_illumina.join(readFileList_nanopore)}
else if (params.filter == '' & params.qual == true){
    readsFileList_qual = readFileList_illumina.join(readFileList_nanopore)}
else if (params.filter == true & params.qual == ''){
    readsFileList_filter = readFileList_illumina.join(readFileList_nanopore)}
else if (params.filter == true & params.qual == true){
    readsFileList_filter = readFileList_illumina.join(readFileList_nanopore)}


readsFileList2_ch.into {readsFileList_ch; test}
test.subscribe{println "got: ${it}"}


process assembly_p {

    input:
        set val(name), file(illumina), file(nanopore) from readsFileList_ch  
    output:
        set val(name), file(illumina), file(nanopore), file('assembly.fasta') into assembly_out2

    script:
    if(params.assembler == 'metaspades')
        """
        spades.py -1 ${illumina[0]} -2 ${illumina[1]}  --meta --nanopore ${nanopore} -o spades_output -t ${task.cpus} -m ${task.memory}
        mv spades_output/contigs.fasta  assembly.fasta

        """
    else if(params.assembler == 'metaflye')
        """
        flye --nano-corr ${nanopore} -o flye_output -t ${task.cpus} --plasmids --meta
        mv flye_output/assembly.fasta assembly.fasta

        """
        //NEED TO ADD THE POLISHING STEP (one round Medaka then Pilon or racon)

}  



assembly_out2.into{assembly_out; test2}
test2.subscribe{println "banana: ${it}"}

if (params.polish=="skip"){assembly_out.set{mapping_ready}}
else {
    process polishing {
        when:
        params.assembler == 'metaflye'

        input:
            set val(name), file(illumina), file(nanopore), file(assembly) from assembly_out
        
        output:
            set val(name), file(illumina), file(nanopore), file('assembly_polished.fasta') into mapping_ready

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
}

//mapping_ready.subscribe{println "skip well done ${it}"}

process assembly_mapping {
    input:
        set val(name), file(illumina), file(nanopore), file(assembly) from mapping_ready

    output:
        set val(name), file(illumina), file(nanopore), file(assembly), file("mapping_illumina.sam"), file("mapping_nanopore.sam") into binning_assembly

    script:
        """
        minimap2 -ax map-ont ${assembly} ${nanopore} > mapping_nanopore.sam
        bwa index -p illumina -a bwtsw ${assembly}
        bwa mem illumina ${illumina} -t ${task.cpus} > mapping_illumina.samm
        
        """
}

//binning_assembly.subscribe{println "mapping done ${it}"}

process binning {
    input:
    set val(name), file(illumina), file(nanopore), file(assembly), file(mapping_illumina), file(mapping_nanopore) from binning_assembly   
        
    output:
        

    script:
    """
    metabat2  -o output_bin -t ${task.cpus} --unbinned --seed 133742 -i ${assembly} ${mapping_illumina} ${mapping_nanopore}

    """ // NEED TO FIND HOW TO ADD MAPPING AND IF RELEVANT
}