#!/usr/bin/env nextflow
nextflow.preview.dsl=2

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

if (params.rna) {
    if (params.dammit_db) {dammit_db=Channel
        .fromPath( params.dammit_db, checkIfExists: true )}
    else {
        include 'modules/dammit_get_databases' params (busco_db : params.busco_db)
        dammit_download_db()
        dammit_db = dammit_download_db.out
    }
}

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

    include 'modules/eggnog'params(output : params.output)
    eggnog(bins_input_ch,eggnog_db)
    bin_annotated_ch=eggnog.out[0].collect().view()
    bin_annot_ch=eggnog.out[0].view()

//************************
// RNA annotation workflow
//************************

// QC
    include fastp_rna from 'modules/fastp'params(output : params.output)
    fastp_rna(rna_input_ch)
    rna_input_ch = fastp_rna.out

// De novo transcript
    include de_novo_transcript from 'modules/trinity_and_salmon'params(output : params.output)
    de_novo_transcript(rna_input_ch)
    transcript_ch=de_novo_transcript.out
// annotation of the transcript
    include 'modules/dammit' params( dammit_user_db : params.dammit_user_db, busco_db: params.busco_db, output: params.output )
    dammit(transcript_ch,dammit_db)
    rna_annotation_ch = dammit.out
// quantification of the annotated transcript
    include quantification from 'modules/trinity_and_salmon'params(output : params.output)
    quantification(rna_input_ch,rna_annotation_ch)
    quant_of_transcrip_ch=quantification.out


//******************************************************
// Parsing bin annot and RNA out into nice graphical out
//******************************************************


// Share pathway to put and HTML file with