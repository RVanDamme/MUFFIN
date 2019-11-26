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

    include 'modules/checkm'params(output : params.output)
    checkm(classify_ch)

//sourmash classification using gtdb database

    include sourmash_bins from 'modules/sourmash'params(output : params.output)
    sourmash_bins(classify_ch,database_sourmash)