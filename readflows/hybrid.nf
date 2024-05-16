#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { modules_to_include } from '../modules/modules_inclusion.nf'
if (params.modular=="full" | params.modular=="assemble" | params.modular=="assem-class" | params.modular=="assem-annot") {
    include {chopper} from '../modules/ont_qc' params(short_qc : params.short_qc, output: params.output)
    include {fastp} from '../modules/fastp' params(output: params.output)
    include {spades} from '../modules/spades' params(output: params.output)
    include {spades_short} from '../modules/spades' params(output: params.output)
    include {flye} from '../modules/flye' params(output : params.output)
    include {minimap_polish} from'../modules/minimap2'
    include {racon} from '../modules/polish'
    include {medaka} from '../modules/polish' params(model : params.model)
    include {pilon} from '../modules/polish' params(output : params.output)
    include {pilon2} from '../modules/polish' params(output : params.output)
    //include {pilong} from '../modules/polish' params(output : params.output)
    include {minimap2} from '../modules/minimap2' //mapping for the binning 
    include {extra_minimap2} from '../modules/minimap2'
    include {bwa} from '../modules/bwa' //mapping for the binning
    include {bwa_bin} from '../modules/bwa' //mapping for the binning
    include {extra_bwa} from '../modules/bwa'   
    include {metabat2} from '../modules/metabat2' params(output : params.output)
    include {metabat2_extra} from '../modules/metabat2' params(output : params.output)
    include {semibin2} from '../modules/semibin2' params(output : params.output, bining_model : params.bining_model, environment : params.environment)
    include {comebin} from '../modules/comebin' params(output : params.output)
    include {separateBins} from '../modules/bins_tools' params(output : params.output)
    include {bin_filter} from '../modules/bins_tools' params(output : params.output)
    include {get_wrong_bin} from '../modules/bins_tools' params(output : params.output)
    include {bin_merger} from '../modules/bins_tools' params(output : params.output)
    include {bam_merger} from '../modules/samtools_merger' params(output : params.output)
    include {bam_merger_extra} from '../modules/samtools_merger' params(output : params.output)
    include {cat_all_bins} from '../modules/cat_all_bins' 
    include {minimap2_bin} from '../modules/minimap2'
    include {metaquast} from '../modules/quast' params(output : params.output)
    include {unmapped_illumina_retrieve} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    include {unmapped_ont_retrieve} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    // include {illumina_reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    // include {ont_reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
    include {sourmash_download_db} from '../modules/sourmashgetdatabase'
    include {sourmash_download_db_full} from '../modules/sourmashgetdatabase'
    include {sourmash_download_db_lineage} from '../modules/sourmashgetdatabase'
    include {checkm_download_db} from '../modules/checkmgetdatabases'
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {
    include {checkm2} from '../modules/checkm2' params(output: params.output, checkm2_low: params.checkm2_low, checkm2db: params.checkm2db)
    include {sourmash_bins} from '../modules/sourmash' params(output: params.output, bintool : params.bintool)
    include {sourmash_ont} from '../modules/sourmash' params(output: params.output)
    include {sourmash_ill} from '../modules/sourmash' params(output: params.output)
    include {sourmash_checkm_parser} from '../modules/checkm_sourmash_parser' params(output: params.output)
    include {get_fasta_path} from '../modules/bins_tools' params(output : params.output)
}
if (params.modular=="full" | params.modular=="annotate" | params.modular=="assem-annot" | params.modular=="class-annot") {
    include {eggnog_download_db} from '../modules/eggnog_get_databases'
    include {eggnog_bin} from '../modules/eggnog' params(output: params.output)
    include {fastp_rna} from '../modules/fastp' params(output: params.output)
    include {de_novo_transcript_and_quant} from '../modules/trinity_and_salmon' params(output: params.output)
    include {eggnog_rna} from '../modules/eggnog' params(output: params.output)
    include {parser_bin_RNA} from '../modules/parser' params(output: params.output)
    include {parser_bin} from '../modules/parser' params(output: params.output)
}

include {readme_output} from '../modules/readme_output' params(output: params.output)
include {test} from '../modules/test_data_dll'

workflow hybrid_workflow{
    if (params.modular=="full" | params.modular=="assemble" | params.modular=="assem-class" | params.modular=="assem-annot") {
        //Channel illumina_input_ch, ont_input_ch

        if (!params.ont || !params.illumina) error "Both ONT and Illumina reads paths must be specified for 'hybride' read type."

        //ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
        ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map {file -> tuple(file.simpleName, file) }.view()
        
        illumina_input_ch = Channel.fromFilePairs("${params.illumina}/*_R{1,2}.fastq{,.gz}", checkIfExists: true)

        if (!params.skip_ont_qc) {
            ont_input_ch = chopper(ont_input_ch)
        }
        if (!params.skip_ill_qc) {
            illumina_input_ch = fastp(illumina_input_ch)
        }

        switch (params.assembler) {
            case "metaspades":
                spades_ch = illumina_input_ch.join(ont_input_ch)
                spades(spades_ch)
                spades.out.set { assembly_ch }
                break

            case "metaflye":
                // MetaFlye suivi d'un polissage hybride
                flye(ont_input_ch)
                // Chaîne de polissage simplifiée
                minimap_polish_ch = minimap_polish(flye.out.join(ont_input_ch))
                racon_ch = racon(ont_input_ch.join(flye.out).join(minimap_polish_ch))
                assembly_ch = medaka(racon_ch)
                //assembly_ch = medaka_ch
                //assembly_ch = pilon(medaka_ch.join(illumina_input_ch), params.polish_iteration)
                break

            default:
                error "Unrecognized assembler: ${params.assembler}. Should be 'metaspades' or 'metaflye'."
        }


        //*********
        // Mapping
        //*********

        // Mapping with Minimap2 for ONT reads
        ont_bam_ch = minimap2(assembly_ch.join(ont_input_ch))
        // Mapping with BWA for Illumina reads
        illumina_bam_ch = bwa(assembly_ch.join(illumina_input_ch))


        // Mapping additional ONT reads if specified
        if (params.extra_ont) {
             extra_ont_ch=Channel.fromPath(params.extra_ont).splitCsv().map { row ->
                        def path = file("${row[0]}")
                        return path
                    }

            extra_bam_ch = extra_minimap2(assembly_ch.join(extra_ont_ch))
        }

        // Mapping additional Illumina reads if specified
        if (params.extra_ill) {
            extra_ill_ch=Channel.fromPath(params.extra_ill).splitCsv().map { row ->
                        def path = file("${row[0]}")
                        return path
                    }

            illumina_extra_bam_ch = extra_bwa(assembly_ch.join(extra_ill_ch))

            if (params.extra_ont){
                extra_bam_ch = bam_merger_extra(extra_bam_ch.join(illumina_extra_bam_ch))
            }
            else{
                extra_bam_ch = illumina_extra_bam_ch
            }
            // Simplified merging logic by using conditional operator
            //extra_bam_ch = params.extra_ont ? bam_merger_extra(extra_bam_ch.join(illumina_extra_bam_ch)) : illumina_extra_bam_ch
        }


        //***************************************************
        // Assembly quality control (metaquast)
        //***************************************************
        if (params.quast) {
            metaquast_out_ch = metaquast(assembly_ch)
        }

        //***************************************************
        // Binning
        //***************************************************
        if (params.bintool) {
        println "Binning tool selected: ${params.bintool}"
        } else {
            println "No binning tool specified, using the default tool: metabat2"
        }


        switch (params.bintool) {
            case 'metabat2':
                if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                    metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                    metabat2_extra(metabat2_ch, extra_bam_ch)
                    bin_out_ch = metabat2_extra.out
                }
                else {
                    metabat2(assembly_ch.join(ont_bam_ch).join(illumina_bam_ch))
                    bin_out_ch = metabat2.out
                }
                break

            case 'semibin2':
                // Prepare bam files
                bam_merger_ch = ont_bam_ch.join(illumina_bam_ch)
                merged_bam_out = bam_merger(bam_merger_ch)
                semibin2_ch = assembly_ch.join(merged_bam_out)
                semibin2(semibin2_ch)
                bin_out_ch = semibin2.out
                break

            case 'comebin':
                // Prepare bam files
                bam_merger_ch = ont_bam_ch.join(illumina_bam_ch)
                merged_bam_out = bam_merger(bam_merger_ch)
                comebin_ch = assembly_ch.join(merged_bam_out)
                comebin(comebin_ch)
                bin_out_ch = comebin.out
                break

            default:
                println "L'outil spécifié (${params.bintool}) n'est pas reconnu. Utilisation de l'outil par défaut: metabat2"
                if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                    metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                    metabat2_extra(metabat2_ch, extra_bam_ch)
                    bin_out_ch = metabat2_extra.out
                }
                else {
                    metabat2(assembly_ch.join(ont_bam_ch).join(illumina_bam_ch))
                    bin_out_ch = metabat2.out
                }
        }
    }
    

    //end of step 1

    //*************************************************
    // STEP 2 classify taxa
    //*************************************************
    if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {

        //**************
        // File handling
        //**************

        //bins (list with id (run id not bin) coma path/to/file)
        if (params.bin_classify) { 
            classify_ch = Channel
                .fromPath( params.bin_classify, checkIfExists: true )
                .splitCsv()
                .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)]  }
                .view()
                }

        else{
            classify_ch = bin_out_ch
            //possibilité de casse !!!
        }


        // Sourmash database configuration
        if (params.sourmash_db) { database_sourmash = file(params.sourmash_db) }
        else {
            sourmash_download_db() 
            database_sourmash = sourmash_download_db.out
        }   
        if (!params.checkm2db){
            println "Checkm2 database configuration"
            // Check that folder and file exist
            path_exists = file(params.db_path).exists() && file("${params.db_path}/${params.db_file}").exists()
            // If the path does not exist or if the file is not found.
            if( !path_exists | params.checkm2db_force_update) {
                checkm_download_db()
            } else {
                println "The specified folder and file already exist : continue."
            }
        }

    
        //*************************
        // Bins classify workflow
        //*************************

        //checkm of the final assemblies
        checkm2(classify_ch, checkm_download_db.out)

        if (!params.skip_bin_sorting){
            bin_filter(checkm2.out.join(classify_ch))
            classify_ch = bin_filter.out
            get_wrong_bin(checkm2.out.join(classify_ch))

            // separateBins(checkm2.out.join(classify_ch))
            // classify_ch = separateBins.out[0]
            // bad_bins_ch = separateBins.out[1]
            
            if (!params.skip_bad_reads_recovery && !params.bin_classify){
                merged_bin_ch = bin_merger(classify_ch)

                bwa_bin_ch = bwa_bin(merged_bin_ch.join(illumina_input_ch))
                minimap_bin_ch = minimap2_bin(merged_bin_ch.join(ont_input_ch))

                // ont_input_ch = ont_reads_retrieval(minimap_bin_ch.join(ont_input_ch))
                // illumina_input_ch = illumina_reads_retrieval(bwa_bin_ch.join(illumina_input_ch))

                unmapped_illumina_retrieve(bwa_bin_ch.join(illumina_input_ch))
                unmapped_ont_retrieve(minimap_bin_ch.join(ont_input_ch))
                if (params.sourmash_db_full) { database_sourmash_full = file(params.sourmash_db_full) }
                else {
                    sourmash_download_db_full() 
                    database_sourmash_full = sourmash_download_db_full.out
                }
                if (params.sourmash_db_lineage) { database_sourmash_lineage = file(params.sourmash_db_lineage) }
                else {
                    sourmash_download_db_lineage() 
                    database_sourmash_lineage = sourmash_download_db_lineage.out
                }
                sourmash_ont(unmapped_ont_retrieve.out, database_sourmash_lineage, database_sourmash_full)
                sourmash_ill(unmapped_illumina_retrieve.out[0], unmapped_illumina_retrieve.out[1], database_sourmash_lineage, database_sourmash_full)
            }
        }
        else {
            classify_ch = get_fasta_path(classify_ch)
        }
         
        if (!params.skip_pilon && params.assembler == 'metaflye' && !params.bin_classify){
            pilon2(classify_ch, illumina_input_ch, params.polish_iteration)
            classify_ch = pilon2.out
        }

        //split the bins for sourmash
        classify_ch.flatMap { name, paths ->
            paths.collect { path -> tuple(name, path) }
        }
        .set { bins_splitted_ch }

        //sourmash classification using gtdb database
        sourmash_bins(bins_splitted_ch,database_sourmash)
        // sourmash_parser_ch = sourmash_bins.out.collect()
    }

    //part 3 
    //*************************************************
    // STEP 3 annotation; kegg pathways + use or RNAseq
    //*************************************************

    if (params.modular=="full" | params.modular=="annotate" | params.modular=="assem-annot" | params.modular=="class-annot") {
        //**************
        // File handling
        //**************

        //RNAseq
        if (params.rna) {rna_input_ch = Channel
                .fromFilePairs( "${params.rna}/*_R{1,2}.fastq{,.gz}", checkIfExists: true)
                .view()
        }

        //bins (list with id (run id not bin) coma path/to/file)
        if (params.bin_annotate) {
            bins_input_ch = Channel
                .fromPath( params.bin_annotate, checkIfExists: true )
                .splitCsv()
                .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)]  }
                .view() 
                }
        else if (params.bin_classify) {
            bins_input_ch = Channel
                .fromPath( params.bin_classify, checkIfExists: true )
                .splitCsv()
                .map { row -> ["${row[0]}", file("${row[2]}", checkIfExists: true)]  }
                .view() 
                }
        else {
            if (params.modular=="assem-annot") {
                println "WARNING : you are using no quality verified bins, be prudent with your results! "
                bins_input_ch = get_fasta_path(bin_out_ch) 
                
            }
            else if (params.modular=="class-annot" || params.modular=="full") {
                bins_input_ch = classify_ch
            }
            else {
                error "No bins given for annotation part."
            }

        }

        //************************
        // Databases Dll and setup
        //************************
        if (params.eggnog_db) {eggnog_db=Channel
                .fromPath( params.eggnog_db, checkIfExists: true )}
        else {
            eggnog_download_db()
            eggnog_db = eggnog_download_db.out
            } 

        //*************************
        // Bins annotation workflow
        //*************************
        bins_input_ch.flatMap { name, paths ->
            paths.collect { path -> tuple(name, path) }
        }
        .set { bins_input_ready_ch }

        eggnog_bin_ch = bins_input_ready_ch.combine(eggnog_db)
        eggnog_bin(eggnog_bin_ch) //annotate the bins
        // bin_annotated_ch=eggnog_bin.out[0].groupTuple(by:0).view()

        //************************
        // RNA annotation workflow
        //************************
        if (params.rna) {
        // QC   
            rna_input_ch = fastp_rna(rna_input_ch) //qc illumina RNA-seq

        // De novo transcript
            de_novo_transcript_and_quant(rna_input_ch) // de novo transcrip assembly and quantification with trinity and salmon
            transcript_ch=de_novo_transcript_and_quant.out
        // annotations of transcript
            eggnog_rna_ch= transcript_ch.combine(eggnog_db)
            eggnog_rna(eggnog_rna_ch) //annotate the RNA-seq transcripts
            rna_annot_ch=eggnog_rna.out[0].view()
        }

        //******************************************************
        // Parsing bin annot and RNA out into nice graphical out
        //******************************************************

        // if (params.rna)  {
        //     parser_bin_RNA(rna_annot_ch,bin_annotated_ch) // parse the annotations in html summary files
        // }
        // else {
        //     parser_bin(bin_annotated_ch) // parse the bins annotation in html summary files
        // }
        // Share pathway to put and HTML file with
    } // end of step 3
    
    readme_output()

   

    

}

workflow.onComplete { 
    log.info ( workflow.success ? "\nDone! Results are stored here --> $params.output \n The Readme file in $params.output describe the structure of the results directories. \n" : "Oops .. something went wrong" )  
}
//***********
//DONE
//***********