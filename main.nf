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
        .fromFilePairs( "${params.illumina}*_R{1,2}.fastq{,.gz}", checkIfExists: true)
        .view() 
// reads_illumina = "${params.illumina}*_R{1,2}.fastq.gz"
// illumina_input_ch = Channel.fromFilePairs(reads_illumina).ifEmpty { error "Cannot find any Illumina reads in the directory: ${params.illumina} \n Delfault is ./illumina \n ${reads_illumina}" }.view()

// reads_ont= "${params.ont}*.fastq.gz"
ont_input_ch = Channel.fromPath("${params.ont}*.fastq{,.gz}",checkIfExists: true).map {file -> tuple(file.simpleName, file) }.view()

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
    --ont                       path to the directory containing the nanopore read file (fastq) (default: ./nanopore)
    --illumina                 path to the directory containing the illumina read file (fastq) (default: ./illumina)
    --rna                       path to the directory containing the illumina read file (fastq) (default: none)

        Output (default output is reassemblies from each bins):
    --output                    path to the output directory (default: $params.output)
    --assembly                  output the original assembly contigs file (default: false)
    --out_qc                    output the reads file after qc (default: false)
    --out_metabat               output the bins produce by metabat2 (default: false)
    --out_concoct               output the bins produce by concoct (default: false)
    --out_maxbin                output the bins produce by meaxbin2 (default: false)
    --out_metawrap              output the bins produce by metawrap refining (default: false)
    --out_bin_reads             output fastq files containing the reads mapped to each bin (default: false)
    --out_unmapped              output sorted bam files containing the unmmaped reads of illumina and nanopore (default:false)


    

        Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB for --metamaps [default: $params.memory]

        Databases:
    --checkm_db                 path to an already INSTALLED checkm database (not the tar file)
    --checkm_tar_db             path to the tar checkm database (it will extract it in the dir)
    --sourmash_db               path to an already installed sourmash database
    --dammit_db                 path to an already installed dammit databases
    --dammmit_user_db           path to a personnal protein database
    --busco_db                  the busco database you want to use in dammit (default: metazoa)
    
        Options:
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
if (params.sourmash_db) { database_sourmash = file(params.sourmash_db) }

else {
    include 'modules/sourmashgetdatabase'
    sourmash_download_db() 
    database_sourmash = sourmash_download_db.out 
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
    include fastp from 'modules/fastp' params(out_qc : params.out_qc, output : params.output)
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
    include sourmash_genome_size from 'modules/sourmash'
    include 'modules/flye' params(assembly : params.assembly, output : params.output)
    include minimap_polish from'modules/minimap2'
    include racon from 'modules/polish'
    include medaka from 'modules/polish' params(model : params.model)
    include pilon from 'modules/polish' params(assembly : params.assembly, output : params.output)
    // FLYE + Pilon 
    flye(sourmash_genome_size(ont_input_ch,database_sourmash))
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
        final_bin_ch = refine2.out
    }
}

else if (params.skip_maxbin2==true) {
    if (  params.skip_metabat2==true || params.skip_concoct==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(concoct_out)
        refine2(refine2_ch, checkm_db_path)
        final_bin_ch = refine2.out
    }
}

else if (params.skip_concoct==true) {
    if (  params.skip_metabat2==true || params.skip_maxbin2==true) {}
    else {
        include refine2 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
        refine2_ch = metabat2_out.join(maxbin2_out)
        refine2(refine2_ch, checkm_db_path)
        final_bin_ch = refine2.out
    }
}

else {
    include refine3 from 'modules/metawrap_refine_bin' params(out_metawrap : params.out_metawrap, output : params.output)
    refine3_ch = metabat2_out.join(maxbin2_out).join(concoct_out)
    refine3(refine3_ch, checkm_db_path)
    final_bin_dir_ch = refine3.out[0]
}

//**************
//Retrieve reads for each bin and assemble them
//**************
if (params.reassembly) {
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
    include unmapped_retrieve from 'modules/seqtk_retrieve_reads'params(output : params.output)
    include 'modules/unicycler_reassemble_from_bin' params(output : params.output)
    include pilon_final from 'modules/polish' params( output : params.output)
    retrieve_unmapped_ch = ill_map_all_bin.join( ont_map_all_bin).join(illumina_input_ch).join(ont_input_ch)
    if (params.out_unmapped == true) {unmapped_retrieve(retrieve_unmapped_ch)}
    retrieve_reads_ch = extract_reads_ch.transpose().combine(ill_map_all_bin, by:0).combine( ont_map_all_bin, by:0).combine(illumina_input_ch, by:0).combine(ont_input_ch, by:0)
    reads_retrieval(retrieve_reads_ch).view()
    unicycler(reads_retrieval.out)

    collected_final_bins_ch=unicycler.out[0].collect()
    final_bins_ch=unicycler.out[0]
}
else {
    collected_final_bins_ch=final_bin_dir_ch
}
//******
// Done
//******



//*************************************************
// STEP 2 classify taxa
//*************************************************

//**************
// File handling
//**************

//bins (list with id (run id not bin) coma path/to/file)
if (params.bin_classify) { 
    classify_ch = Channel
        .fromPath( params.bin_classify, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}","${row[1]}", file("${row[2]}", checkIfExists: true)]  }
        .view()
        }

else {classify_ch=final_bins_ch}

//*************************
// Bins classify workflow
//*************************

//checkm of the final assemblies

    include 'modules/checkm'params(output : params.output)
    checkm(classify_ch.groupTuple(by:0))

//sourmash classification using gtdb database

    include sourmash_bins from 'modules/sourmash'params(output : params.output)
    sourmash_bins(classify_ch,database_sourmash)

    include sourmash_checkm_parser from 'modules/checkm_sourmash_parser'params(output: params.output)
    sourmash_checkm_parser(checkm.out[0],sourmash_bins.out.collect())


//*************************************************
// STEP 3 annotation; kegg pathways + use or RNAseq
//*************************************************

//**************
// File handling
//**************

//RNAseq
if (params.rna) {rna_input_ch = Channel
        .fromFilePairs( "${params.rna}*_R{1,2}.fastq{,.gz}", checkIfExists: true)
        .view()
}

//bins (list with id (run id not bin) coma path/to/file)
if (params.bin_annotate) {
    bins_input_ch = Channel
        .fromPath( params.bin_annotate, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}","${row[1]}", file("${row[2]}", checkIfExists: true)]  }
        .view() 
        }
else {bins_input_ch = final_bins_ch }



//************************
// Databases Dll and setup
//************************
if (params.eggnog_db) {eggnog_db=Channel
        .fromPath( params.eggnog_db, checkIfExists: true )}
else {
    include 'modules/eggnog_get_databases'
    eggnog_download_db()
    eggnog_db = eggnog_download_db.out
    }
//*************************
// Bins annotation workflow
//*************************

    include eggnog_bin from 'modules/eggnog'params(output : params.output)
    eggnog_bin_ch= bins_input_ch.combine(eggnog_db)
    eggnog_bin(eggnog_bin_ch)
    bin_annotated_ch=eggnog_bin.out[0].groupTuple(by:0).view()

//************************
// RNA annotation workflow
//************************

// QC
    include fastp_rna from 'modules/fastp'params(output : params.output)
    fastp_rna(rna_input_ch)
    rna_input_ch = fastp_rna.out

// De novo transcript
    include de_novo_transcript_and_quant from 'modules/trinity_and_salmon'params(output : params.output)
    de_novo_transcript_and_quant(rna_input_ch)
    transcript_ch=de_novo_transcript_and_quant.out
// annotations of transcript
    include eggnog_rna from 'modules/eggnog'params(output : params.output)
    eggnog_rna_ch= transcript_ch.combine(eggnog_db)
    eggnog_rna(eggnog_rna_ch)
    rna_annot_ch=eggnog_rna.out[0].view()


//******************************************************
// Parsing bin annot and RNA out into nice graphical out
//******************************************************

    include from 'modules/parser'params(output: params.output)
    parser(rna_annot_ch,bin_annotated_ch)

// Share pathway to put and HTML file with

//***********
// MAFIN DONE
//***********