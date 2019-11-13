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
if (params.list_bins) {bins_pre_input_ch = Channel
        .fromPath( params.list_bins, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", "${row[1]}", file("${row[2]}", checkIfExists: true)]  }
        .view() 
        
        include fasta_check from 'modules/ext_check'
        fasta_check(bins_pre_input_ch)
        bins_input_ch=fasta_check.out
        }

else if () {bins_input_ch = unicycler.out[0] }



//************************
// Databases Dll and setup
//************************

if (params.rna) {
    if (params.dammmit_db) {dammit_db=params.dammmit_db}
    else {
        include 'modules/dammit_get_databases'
        dammit_download_db()
        dammit_db = dammit_download_db.out
    }
}

if (params.eggnog_db) {eggnog_db=params.eggnog_db}
else {
    include 'modules/eggnog_get_databases'
    eggnog_download_db()
    eggnog_db = eggnog_download_db.out
    }

//*************************
// Bins annotation workflow
//*************************

include 'modules/eggnog'
eggnog(bins_input_ch,eggnog_db)
bin_annotated_ch=eggnog.out

//************************
// RNA annotation workflow
//************************

// QC
include fastp_rna from 'modules/fastp'
    fastp(rna_input_ch)
    rna_input_ch = fastp.out

// De novo transcript
include from 'modules/trinity_and_salmon'
    de_novo_transcript_and_quant(rna_input_ch)
    transcript_ch=de_novo_transcript_and_quant.out[0]
    quant_of_transcrip_ch=de_novo_transcript_and_quant.out[1]
// annotation of the transcript

include from 'modules/dammit' params(output : params.output, dammmit_user_db : params.dammmit_user_db )
    dammit(transcript_ch,dammit_db)
    rna_annotation_ch = dammit.out

//******************************************************
// Parsing bin annot and RNA out into nice graphical out
//******************************************************


// Share pathway to put and HTML file with