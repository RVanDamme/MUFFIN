#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { modules_to_include } from '../modules/modules_inclusion.nf'
if (params.modular=="full" | params.modular=="assemble" | params.modular=="assem-class" | params.modular=="assem-annot") {
    include {chopper} from '../modules/ont_qc' params(output: params.short_qc)
    include {fastp} from '../modules/fastp' params(output: params.output)
    include {spades} from '../modules/spades' params(output: params.output)
    include {spades_short} from '../modules/spades' params(output: params.output)
    include {flye} from '../modules/flye' params(output : params.output)
    include {minimap_polish} from'../modules/minimap2'
    include {racon} from '../modules/polish'
    include {medaka} from '../modules/polish' params(model : params.model)
    include {pilon} from '../modules/polish' params(output : params.output)
    include {pilong} from '../modules/polish' params(output : params.output)
    include {minimap2} from '../modules/minimap2' //mapping for the binning 
    include {extra_minimap2} from '../modules/minimap2'
    include {bwa} from '../modules/bwa' //mapping for the binning
    include {extra_bwa} from '../modules/bwa'
    //include {metabat2_extra} from '../modules/metabat2' params(output : params.output)    
    include {metabat2} from '../modules/metabat2' params(output : params.output)
    include {semibin2} from '../modules/semibin2' params(output : params.output)
    include {comebin} from '../modules/comebin' params(output : params.output)
    include {bam_merger} from '../modules/samtools_merger' params(output : params.output)
    //include {contig_list} from '../modules/list_ids'
    include {cat_all_bins} from '../modules/cat_all_bins'
    include {bwa_bin} from '../modules/bwa'  
    include {minimap2_bin} from '../modules/minimap2'
    include {metaquast} from '../modules/quast' params(output : params.output)
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
    include {sourmash_download_db} from '../modules/sourmashgetdatabase'
    include {checkm_download_db} from '../modules/checkmgetdatabases'
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {
    include {reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    include {checkm2} from '../modules/checkm2' params(output: params.output)
    include {sourmash_bins} from '../modules/sourmash' params(output: params.output)
    include {sourmash_checkm_parser} from '../modules/checkm_sourmash_parser' params(output: params.output)
    //include {unmapped_retrieve} from '../modules/seqtk_retrieve_reads' params(output : params.output)
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
params.db_path = '/home/user/databases'
params.db_file = 'uniref100.KO.1.dmnd'


workflow long_read_workflow{
    // Initialisation des variables pour les chemins des bases de données
    if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){

        // Configuration de la base de données Sourmash
        if (params.sourmash_db) { database_sourmash = file(params.sourmash_db) }
        else {
            sourmash_download_db() 
            database_sourmash = sourmash_download_db.out
        }   
        if (!params.checkm2db){
            // Vérification de l'existence du dossier et du fichier
            path_exists = file(params.db_path).exists() && file("${params.db_path}/${params.db_file}").exists()

            // Si le chemin n'existe pas ou si le fichier n'est pas trouvé.
            if( !path_exists | params.checkm2db_force_update) {
                checkm_download_db()
            } else {
                println "Le dossier et le fichier spécifié existent déjà."
            }
        }
    }

    // Assemblage avec des reads longs (ONT)
    if (!params.ont) error "ONT reads path must be specified for 'long' read type."
    ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
    if (!params.skip_ont_qc) {
        ont_input_ch = chopper(ont_input_ch)
    }

    // Assemblage avec des reads longs (ONT) utilisant Flye
    flye(ont_input_ch)
    //assembly_ch = flye.out

    minimap_polish_ch = minimap_polish(flye.out.join(ont_input_ch))
    racon_ch = racon(ont_input_ch.join(flye.out).join(minimap_polish_ch))
    medaka_ch = medaka(racon_ch)
    assembly_ch = pilong(medaka_ch.join(ont_input_ch), params.polish_iteration)

    // assembly_ch.flatMap { contigs ->
    //     minimap_polish(contigs, ont_input_ch)
    //     }.flatMap { polished ->
    //         racon(polished)
    //     }.flatMap { racon_out ->
    //         medaka(racon_out)
    //     }.flatMap { medaka_out ->
    //         pilong(medaka_out.join(ont_input_ch), params.polish_iteration)
    //     }.set { assembly_ch }
    //     //ajouter pilon

    //*********
    // Mapping
    //*********
    
    // Mapping with Minimap2 for ONT reads
    // ont_bam_ch = assembly_ch
    //     .join(ont_input_ch)
    //     .map { assembly, reads -> [assembly, reads] }
    //     .flatMap { minimap2(it) }

    // // Mapping additional ONT reads if specified
    // if (params.extra_ont) {
    //     ont_extra_bam_ch = assembly_ch
    //         .join(Channel.fromPath(params.extra_ont))
    //         .map { assembly, extraReads -> [assembly, extraReads] }
    //         .flatMap { extra_minimap2(it) }
    // }

    // Mapping with Minimap2 for ONT reads
    ont_bam_ch = minimap2(assembly_ch.join(ont_input_ch))

    // Mapping additional ONT reads if specified
    if (params.extra_ont) {
        ont_extra_bam_ch = minimap2(assembly_ch.join(extra_ont_ch))
    }

    //***************************************************
    // Assembly quality control
    //***************************************************
    if (params.reference) {
        //Channel ref_ch = params.reference
        ref_ch = Channel.fromPath(params.reference)
        metaquast(assembly_ch.join(ref_ch))
        metaquast_out_ch = metaquast.out
    }

    //***************************************************
    // Binning
    //***************************************************
    if (params.bintool) {
    println "Outil de binning sélectionné: ${params.bintool}"
    } else {
        println "Aucun outil de binning spécifié, utilisation de l'outil par défaut: metabat2"
    }

    switch (params.bintool) {
        case 'metabat2':
            if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                metabat2_ch = ont_bam_ch.join(extra_bam)
                bam_merger(metabat2_ch)
                metabat2_extra(assembly_ch, bam_merger.out)
                metabat2_out = metabat2_extra.out
            }
            else {
                metabat2_ch = assembly_ch.join(ont_bam_ch)
                metabat2(metabat2_ch)
                metabat2_out = metabat2.out
            }
            break

        case 'semibin2':
            semibin2_ch = assembly_ch.join(ont_bam_ch)
            semibin2(semibin2_ch)
            semibin2_out = semibin2.out
            break

        case 'comebin':
            comebin_ch = assembly_ch.join(ont_bam_ch)
            comebin(comebin_ch)
            comebin_out = comebin.out
            break

        default:
            println "L'outil spécifié (${params.bintool}) n'est pas reconnu. Utilisation de l'outil par défaut: metabat2"
            if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                metabat2_ch = ont_bam_ch.join(extra_bam)
                bam_merger(metabat2_ch)
                metabat2_extra(assembly_ch, bam_merger.out)
                metabat2_out = metabat2_extra.out
            }
            else {
                metabat2_ch = assembly_ch.join(ont_bam_ch)
                metabat2(metabat2_ch)
                metabat2_out = metabat2.out
            }
            
    }


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


        //merge every bining tool result
        else {classify_ch= metabat2_out_ch.join(semibin2_out_ch).join(comebin_out_ch)}

        // if (params.modular=="classify" | params.modular=="class-annot") {
        //     // sourmash_db
        //     if (params.sourmash_db) { database_sourmash = file(params.sourmash_db) }
        //     else {
        //         sourmash_download_db() 
        //         database_sourmash = sourmash_download_db.out
        //     }   
        //     if (!params.checkm2db){
        //         // Vérification de l'existence du dossier et du fichier
        //         path_exists = file(params.db_path).exists() && file("${params.db_path}/${params.db_file}").exists()

        //         // Si le chemin n'existe pas ou si le fichier n'est pas trouvé.
        //         if( !path_exists | params.checkm2db_force_update) {
        //             checkm_download_db()
        //         } else {
        //             println "Le dossier et le fichier spécifié existent déjà."
        //         }
        //     }
        // }
    
        //*************************
        // Bins classify workflow
        //*************************

        //checkm of the final assemblies
        //checkm(classify_ch.groupTuple(by:0)) //checkm QC of the bins
        checkm2(classify_ch)
        checkm2_out_ch = checkm2.out 

        //sourmash classification using gtdb database
        sourmash_bins(classify_ch,database_sourmash) // fast classification using sourmash with the gtdb (not the best classification but really fast and good for primarly result)
        sourmash_checkm_parser(checkm.out[0],sourmash_bins.out.collect()) //parsing the result of sourmash and checkm in a single result file
    }
}



