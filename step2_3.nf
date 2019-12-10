#!/usr/bin/env nextflow
nextflow.preview.dsl=2

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

else {classify_ch=all_bins_ch}



//************************
// Databases Dll and setup
//************************


    include 'modules/sourmashgetdatabase'
    sourmash_download_db() 
    database_sourmash = sourmash_download_db.out 


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


//*************************
// Bins classify workflow
//*************************

//checkm of the final assemblies

    include 'modules/checkm' params (output : params.output)
    checkm(classify_ch.groupTuple(by:0))

//sourmash classification using gtdb database

    include sourmash_bins from 'modules/sourmash'params(output : params.output)
    sourmash_bins(classify_ch,database_sourmash)

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
if (params.bin_classify) {
    bins_input_ch = Channel
        .fromPath( params.bin_classify, checkIfExists: true )
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
    include parser from 'modules/parser'params(output: params.output)
    parser(rna_annot_ch,bin_annotated_ch)

// Share pathway to put and HTML file with

//***********
// MAFIN DONE
//***********
