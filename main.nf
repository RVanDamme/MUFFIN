#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// if params.mode==retrieval {whole step one script}
// if params.mode==analysis {whole bin analysis and genome analysis}
// TODO publish of file in output (need to deicide what to keep)


//*****************************************
// Input, check(command) and stdout (params and such)
//*****************************************

// fool proof checks
if (params.help) { exit 0, helpMSG() }

if (params.assembler!='metaflye' && params.assembler!='metaspades') {
    exit 1, "--assembler: ${params.method}. Should be 'metaflye' or 'metaspades'"}

// stdout early usage (print header + default or modified params)

// DATA INPUT (ONT and ILLUMINA)
illumina_input_ch = Channel
        .fromFilePairs( "${params.illumina}*_R{1,2}.fastq.gz", checkIfExists: true)
        .view() 
// reads_illumina = "${params.illumina}*_R{1,2}.fastq.gz"
// illumina_input_ch = Channel.fromFilePairs(reads_illumina).ifEmpty { error "Cannot find any Illumina reads in the directory: ${params.illumina} \n Delfault is ./illumina \n ${reads_illumina}" }.view()

// reads_ont= "${params.ont}*.fastq.gz"
ont_input_ch = Channel.fromPath("${params.ont}*.fastq.gz",checkIfExists: true).map {file -> tuple(file.simpleName, file) }.view()

/*
// Help Message
def helpMSG() {
    log.info """
    *********Metagenomic Assembly pipeline using nextFlow for Illumina and Nanopore reads*********

    Mafin is composed of 2 part the retrieval of potential genome and the analysis of said genomes

        Usage example for retrieval:
    nextflow run mafin --retrieve --ont /path/to/ont_dir --illumina /path/to/illumina_dir --metaspades -profile conda
    or 
    nextflow run mafin --retrieve --ont /path/to/ont_dir --illumina /path/to/illumina_dir --metaflye -profile conda

        Input:
    --ont                       path to the directory containing the nanopore read file (fastq)
    -- illumina                 path to the directory containing the illumina read file (fastq)

        Output (default output is reassemblies from each bins):
    --output                    path to the output directory (default: $params.output)
    --assembly                  output the original assembly contigs file (default: false)
    --out_qc                    output the reads file after qc (default: false)
    --out_metabat               output the bins produce by metabat2 (default: false)
    --out_concoct               output the bins produce by concoct (default: false)
    --out_maxbin                output the bins produce by meaxbin2 (default: false)
    --out_metawrap              output the bins produce by metawrap refining (default: false)
    --out_bin_reads             output fastq files containing the reads map to each bin (default: false)


    

        Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB for --metamaps [default: $params.memory]
    
        Options:
    --checkm_db                 path to an already INSTALLED checkm database (not the tar file)
    --checkm_tar_db             path to the tar checkm database (it will extract it in the dir)
    --sourmash                  path to an already installed sourmash database
    --skip_ill_qc               skip quality control of illumina files
    --skip_ont_qc               skip quality control of nanopore file
    --short_qc                  minimum size of the reads to be kept (default: $params.short_qc )
    --filtlong                  use filtlong to improve the quality furthermore (default: false)
    --polish_iteration          number of iteration of the polish step (advanced)
    --polish_threshold          threshold to reach to stop the iteration of the polish step (advanced)
    --skip_metabat2             skip the binning using metabat2 (advanced)
    --skip_maxbin2              skip the binning using maxbin2 (advanced)
    --skip_concoct              skip the binning using concoct (advanced)

        Nextflow options:
    -profile                    change the profile of nextflow (currently available conda)
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)
    """
}
/***********************
* Whole process using MODULES
************************/

// Loading modules

//******************************
// Databases check or retrieval
//******************************

// sourmash_db
if (params.assembler=="metaflye") { 
    if (params.sour_db) { database_sourmash = file(params.sour_db) }

    else {
        include 'modules/sourmashgetdatabase'
        sourmash_download_db() 
        database_sourmash = sourmash_download_db.out 
    }
}

// checkm_db
if (params.checkm_db) {
    include 'modules/checkmsetupDB'
    untar = true
    checkm_setup_db(params.checkm_db, untar)
}

else if (params.checkm_tar_db) {
    include 'modules/checkmsetupDB'
    untar = false
    checkm_setup_db(params.checkm_db, untar)
}

else {
    include 'modules/checkmsetupDB'
    include 'modules/checkmgetdatabases'
    untar = false
    checkm_setup_db(checkm_download_db(), untar)
}


if (params.skip_ont_qc == true) {}
else if (params.skip_ont_qc==false){
    include discard_short from 'modules/ont_qc' params(short_qc : params.short_qc)
    split_ont_ch = ont_input_ch.splitFastq(by:100000, file:true)
    discard_short(split_ont_ch)
    if (params.filtlong==true){
        include filtlong from 'modules/ont_qc'
        filtlong(discard_short.out)
        merging_ch = filtlong.out.groupTuple()
    }
    else {
        merging_ch = discard_short.out.groupTuple()
    }
    include merge from 'modules/ont_qc' params(out_qc : params.out_qc, output : params.output)
    merge(merging_ch) 
    ont_input_ch = merge.out

}


// QC check Illumina

if (params.skip_ill_qc==true) {}

else if (params.skip_ill_qc==false){
    include 'modules/fastp' params(out_qc : params.out_qc, output : params.output)
    fastp(illumina_input_ch)
    illumina_input_ch = fastp.out
}

//**********
// Assembly 
//**********

// Meta-SPADES

if (params.assembler=="metaspades") {
    include 'modules/spades' params(assembly : params.assembly, output : params.output)
    spades_ch= illumina_input_ch.join(ont_input_ch)
    spades(spades_ch)
    assembly_ch = spades.out
}


// Meta-FLYE

if (params.assembler=="metaflye") {
    include 'modules/flye' params(assembly : params.assembly, output : params.output)
    include 'modules/pilon' params(assembly : params.assembly, output : params.output)
    include minimap2 as minimap_polish from'modules/minimap2' // determine if it's ONT, ILL or both to map for pilon
    // FLYE + Pilon 
    pilon(minimap2_polish(flye(ont_input_ch), ont_input_ch), flye.out, params.polish_iteration, params.polish_threshold) // don't remember which reads to map
    assembly_ch = pilon.out
}

//*********
// Mapping
//*********

// ONT mapping

    include minimap2 from 'modules/minimap2'
    minimap2_ch = assembly_ch.join(ont_input_ch)
    minimap2(minimap2_ch)
    ont_bam_ch = minimap2.out


// Illumina mapping

    include bwa from 'modules/bwa'
    bwa_ch = assembly_ch.join(illumina_input_ch)
    bwa(bwa_ch)
    illumina_bam_ch = bwa.out

//***************************************************
// Binning
//***************************************************

// metabat2 + checkm  OR metabat2 and check of it separately?

if (params.skip_metabat2==true) {}

else {
    include 'modules/metabat2' params(out_metabat : params.out_metabat, output : params.output)
    metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
    metabat2(metabat2_ch)
    metabat2_out = metabat2.out
}

// Maxbin2 OR CheckM Maxbin2

if (params.skip_maxbin2==true) {}

else {
    include 'modules/maxbin2' params(out_maxbin : params.out_maxbin, output : params.output)
    maxbin2_ch = assembly_ch.join(ont_input_ch).join(illumina_input_ch)
    maxbin2(maxbin2_ch)
    maxbin2_out = maxbin2.out
}

// Concoct OR CheckM Concoct

if (params.skip_concoct==true) {}

else {
    include 'modules/concoct' params(out_concoct : params.out_concoct, output : params.output)
    concoct_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
    concoct(concoct_ch)
    concoct_out = concoct.out
}

// Bin refine

if (params.skip_metabat2==true) {
    if (  params.skip_maxbin2==true || params.skip_concoct==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = maxbin2_out.join(concoct_out)
        refine2(refine2_ch)
        final_bin_ch = refine2.out
    }
}

else if (params.skip_maxbin2==true) {
    if (  params.skip_metabat2==true || params.skip_concoct==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(concoct_out)
        refine2(refine2_ch)
        final_bin_ch = refine2.out
    }
}

else if (params.skip_concoct==true) {
    if (  params.skip_metabat2==true || params.skip_maxbin2==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(maxbin2_out)
        refine2(refine2_ch)
        final_bin_ch = refine2.out
    }
}

else {
    include refine3 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
    refine3_ch = metabat2_out.join(maxbin2_out).join(concoct_out)
    refine3(refine3_ch)
    final_bin_ch = refine3.out
}

//**************
//Retrieve reads for each bin and assemble them
//**************

// retrieve the ids of each bin contigs

    include contig_list from 'modules/list_ids'
    contig_list(final_bin_ch)
    extract_reads_ch = contig_list.out.transpose()

 
// bam align the reads to ALL OF THE CONTIGS 

    include 'modules/cat_all_bins'
    include bwa_bin from 'modules/bwa'  
    include minimap2_bin from 'modules/minimap2'
    cat_all_bins(final_bin_ch)
    fasta_all_bin = cat_all_bins.out
    bwa_all_bin = fasta_all_bin.join(illumina_input_ch)
    ill_map_all_bin = bwa_bin(bwa_all_bin)    
    minimap2_all_bin = fasta_all_bin.join(ont_input_ch)
    ont_map_all_bin = minimap2_bin(minimap2_all_bin)

// retrieve the reads aligned to the contigs + run unicycler

    include reads_retrieval from 'modules/list_ids'params(out_bin_reads: params.out_bin_reads, output : params.output)
    include 'modules/unicycler_reassemble_from_bin' params(output : params.output)
    retrieve_reads_ch = extract_reads_ch.join(ill_map_all_bin).join( ont_map_all_bin).join(illumina_input_ch).join(ont_input_ch)
    unicycler(reads_retrieval(retrieve_reads_ch))
    
//******
// Done
//******
