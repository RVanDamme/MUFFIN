#!/usr/bin/env nextflow
nextflow.preview.dsl=2

// if params.mode==retrieval {whole step one script}
// if params.mode==analysis {whole bin analysis and genome analysis}
// TODO publish of file in output (need to deicide what to keep)


//*****************************************
// Input, check(command) and stdout (params and such)
//*****************************************
// if (params.out_all == true) {
//     params.assembly = true
//     params.out_qc = true
//     params.out_metawrap = true
//     params.out_bin_reads = true
//     params.out_unmapped = true
//     }
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

// extra ont reads
if (params.extra_ont != false) {
extra_ont_ch=Channel.fromPath(params.extra_ont).splitCsv().map { row ->
            def path = file("${row[0]}")
            return path
        }
}

// extra ill reads
if (params.extra_ill != false) {
extra_ill_ch=Channel.fromPath(params.extra_ill).splitCsv().map { row ->
            def path = file("${row[0]}")
            return path
        }
}

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
        Outputed files:
    assembly                    the original assembly contigs file 
    qc                          the reads file after qc
    metabat                     the bins produce by metabat2
    concoct                     the bins produce by concoct
    maxbin                      the bins produce by meaxbin2
    metawrap                    the bins produce by metawrap refining
    mapped bin reads            the fastq files containing the reads mapped to each metawrap bin
    unmapped bin reads          the fastq files containing the unmmaped reads of illumina and nanopore


    

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
    --model                     the model medaka will use (default: r941_min_high)
    --polish_iteration          number of iteration of pilon in the polish step (advanced)
    --extra_ill                 a list of additional ill sample file (with full path with a * instead of _R1,2.fastq) to use for the binning in Metabat2 and concoct
    --extra_ont                 a list of additional ont sample file (with full path) to use for the binning in Metabat2 and concoct
    --SRA_ill                   a list of additional ill sample from SRA accession number to use for the binning in Metabat2 and concoct
    --SRA_ont                   a list of additional ont sample from SRA accession number to use for the binning in Metabat2 and concoct
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

//**********************************
// Databases and metawrap obtention
//**********************************

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
    checkm_db_path = checkm_setup_db.out
}

else if (params.checkm_tar_db) {
    include 'modules/checkmsetupDB'
    untar = false
    checkm_setup_db(params.checkm_db, untar)
    checkm_db_path = checkm_setup_db.out
}

else {
    include 'modules/checkmsetupDB'
    include 'modules/checkmgetdatabases'
    untar = false
    checkm_setup_db(checkm_download_db(), untar)
    checkm_db_path = checkm_setup_db.out
}

//************
// QC OF READS
//************
if (params.skip_ont_qc == true) {}
else if (params.skip_ont_qc==false){
    include discard_short from 'modules/ont_qc' params(short_qc : params.short_qc)
    split_ont_ch = ont_input_ch.splitFastq(by:100000, file:true)
    discard_short(split_ont_ch)
    if (params.filtlong==true){
        include filtlong from 'modules/ont_qc' params(short_qc : params.short_qc)
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
    include 'modules/sourmash'
    include 'modules/flye' params(assembly : params.assembly, output : params.output)
    include minimap_polish from'modules/minimap2'
    include racon from 'modules/polish'
    include medaka from 'modules/polish' params(model : params.model)
    include pilon from 'modules/polish' params(assembly : params.assembly, output : params.output)
    // FLYE + Pilon 
    flye(sourmash(ont_input_ch,database_sourmash))
    flye_to_map = flye.out.join(ont_input_ch)
    minimap_polish(flye_to_map)
    map_to_racon = ont_input_ch.join(flye.out).join(minimap_polish.out)
    medaka(racon(map_to_racon))
    medaka_to_pilon = medaka.out.join(ont_input_ch)
    pilon(medaka_to_pilon, params.polish_iteration)
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

// Illumina mapping of extra files

if (params.extra_ont != false) {
    include extra_minimap2 from 'modules/minimap2'
    minimap_extra = assembly_ch.join(extra_ont_ch)
    extra_minimap2(minimap_extra)
    ont_extra_bam = extra_minimap2.out.collect()
}

// Illumina mapping

    include bwa from 'modules/bwa'
    bwa_ch = assembly_ch.join(illumina_input_ch)
    bwa(bwa_ch)
    illumina_bam_ch = bwa.out

// Illumina mapping of extra files

if (params.extra_ill != false) {
    include extra_bwa from 'modules/bwa'
    bwa_extra = assembly_ch.join(extra_ill_ch)
    extra_bwa(bwa_extra)
    illumina_extra_bam = extra_bwa.out.collect()
}

if (params.extra_ont != false && params.extra_ill != false ) {
    extra_bam = illumina_extra_bam.concat(ont_extra_bam)
}
else if (params.extra_ont != false) {
    extra_bam = ont_extra_bam
}
else if (params.extra_ill != false) {
    extra_bam = illumina_extra_bam
}
//***************************************************
// Binning
//***************************************************

// metabat2 + checkm  OR metabat2 and check of it separately?

if (params.skip_metabat2==true) {}

else {
    if (params.extra_ont != false || params.extra_ill != false ) {
        include metabat2_extra from 'modules/metabat2' params(out_metabat : params.out_metabat, output : params.output)
        metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        metabat2_extra(metabat2_ch, extra_bam)
        metabat2_out = metabat2_extra.out
    }
    else {    
        include metabat2 from 'modules/metabat2' params(out_metabat : params.out_metabat, output : params.output)
        metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        metabat2(metabat2_ch)
        metabat2_out = metabat2.out
    }
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
    if (params.extra_ont != false || params.extra_ill != false ) {
        include concoct_extra from 'modules/concoct' params(out_concoct : params.out_concoct, output : params.output)
        concoct_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        concoct_extra(concoct_ch, extra_bam)
        concoct_out = concoct_extra.out
    }
    else {
        include concoct from 'modules/concoct' params(out_concoct : params.out_concoct, output : params.output)
        concoct_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        concoct(concoct_ch)
        concoct_out = concoct.out
    }
}

// Bin refine

if (params.skip_metabat2==true) {
    if (  params.skip_maxbin2==true || params.skip_concoct==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = maxbin2_out.join(concoct_out)
        refine2(refine2_ch, checkm_db_path)
        final_bin_ch = refine2.out[0]
    }
}

else if (params.skip_maxbin2==true) {
    if (  params.skip_metabat2==true || params.skip_concoct==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(concoct_out)
        refine2(refine2_ch, checkm_db_path)
        final_bin_ch = refine2.out[0]
    }
}

else if (params.skip_concoct==true) {
    if (  params.skip_metabat2==true || params.skip_maxbin2==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(maxbin2_out)
        refine2(refine2_ch, checkm_db_path)
        final_bin_ch = refine2.out[0]
    }
}

else {
    include refine3 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
    refine3_ch = metabat2_out.join(maxbin2_out).join(concoct_out)
    refine3(refine3_ch, checkm_db_path)
    final_bin_ch = refine3.out[0]
}

//**************
//Retrieve reads for each bin and assemble them
//**************

// retrieve the ids of each bin contigs

    include 'modules/list_ids'
    contig_list(final_bin_ch)
    extract_reads_ch = contig_list.out.view()

 
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

// retrieve the reads aligned to the contigs + run unicycler + polish with pilon for 2 round

    include reads_retrieval from 'modules/seqtk_retrieve_reads'params(out_bin_reads: params.out_bin_reads, output : params.output)
    include unmapped_retrieve from 'modules/seqtk_retrieve_reads'params(out_unmapped: params.out_unmapped, output : params.output)
    include 'modules/unicycler_reassemble_from_bin' params(output : params.output)
    retrieve_unmapped_ch = ill_map_all_bin.join( ont_map_all_bin).join(illumina_input_ch).join(ont_input_ch)
    if (params.out_unmapped == true) {unmapped_retrieve(retrieve_unmapped_ch)}
    retrieve_reads_ch = extract_reads_ch.join(ill_map_all_bin).join( ont_map_all_bin).join(illumina_input_ch).join(ont_input_ch).transpose()
    reads_retrieval(retrieve_reads_ch).view()
    unicycler(reads_retrieval.out)
    final_assemblies_ch=unicycler.out[0].collect()
//checkm of the final assemblies

    include 'modules/checkm'params(output : params.output)
    checkm(final_assemblies_ch)

//******
// Done
//******



//*******************************************************************************
// STEP 2 Taxonomy; annotation; kegg pathways (and maybe go-term) + use or RNAseq
//*******************************************************************************
