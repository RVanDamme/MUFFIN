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

workflow hybrid_workflow{
    // Initialisation des variables pour les chemins des bases de données
    if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
        Channel database_sourmash

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
    

    //Channel illumina_input_ch, ont_input_ch

    if (!params.ont || !params.illumina) error "Both ONT and Illumina reads paths must be specified for 'hybride' read type."

    //ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
    ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map {file -> tuple(file.simpleName, file) }.view()

    illumina_input_ch = Channel.fromFilePairs("${params.illumina}/*_R{1,2}.fastq{,.gz}", checkIfExists: true)

    if (!params.skip_ont_qc) {
        ont_input_ch = ont_input_ch.flatMap { chopper(it) }
    }
    if (!params.skip_ill_qc) {
        illumina_input_ch = illumina_input_ch.flatMap { fastp(it) }
    }


    //Channel assembly_ch

    switch (params.assembler) {
        case "metaspades":
            spades_ch = illumina_input_ch.join(ont_input_ch)

            // Run spades process with mixed reads
            spades(spades_ch)
            spades.out.set { assembly_ch }
            break

        case "metaflye":
            // MetaFlye suivi d'un polissage hybride
            flye(ont_input_ch)
            assembly_ch = flye.out
            // Chaîne de polissage simplifiée
            assembly_ch.flatMap { contigs ->
                minimap_polish(contigs, ont_input_ch)
            }.flatMap { polished ->
                racon(polished)
            }.flatMap { racon_out ->
                medaka(racon_out)
            }.flatMap { medaka_out ->
                pilon(medaka_out, illumina_input_ch, params.polish_iteration)
            }.set { assembly_ch }
            break

        default:
            error "Unrecognized assembler: ${params.assembler}. Should be 'metaspades' or 'metaflye'."
    }


    //*********
    // Mapping
    //*********
    
    // Mapping with Minimap2 for ONT reads
    ont_bam_ch = assembly_ch
        .join(ont_input_ch)
        .map { assembly, reads -> [assembly, reads] }
        .flatMap { minimap2(it) }

    // Mapping additional ONT reads if specified
    if (params.extra_ont) {
        ont_extra_bam_ch = assembly_ch
            .join(Channel.fromPath(params.extra_ont))
            .map { assembly, extraReads -> [assembly, extraReads] }
            .flatMap { extra_minimap2(it) }
    }

    // Mapping with BWA for Illumina reads
    illumina_bam_ch = assembly_ch
        .join(illumina_input_ch)
        .map { assembly, reads -> [assembly, reads] }
        .flatMap { bwa(it) }

    // Mapping additional Illumina reads if specified
    if (params.extra_ill) {
        illumina_extra_bam_ch = assembly_ch
            .join(Channel.fromPath(params.extra_ill))
            .map { assembly, extraReads -> [assembly, extraReads] }
            .flatMap { extra_bwa(it) }
    }


    //***************************************************
    // Assembly quality control-> may be pushed in a function
    //***************************************************
    if (params.reference) {
        //Channel ref_ch = params.reference
        ref_ch = Channel.fromPath(params.reference)
        // metaquast(assembly_ch, ref_ch)
        // metaquast_out_ch = metaquast.out
        metaquast_out_ch = metaquast(assembly_ch, ref_ch)
    }

    //***************************************************
    // Binning
    //***************************************************
    if (params.bintool) {
    println "Outil de binning sélectionné: ${params.bintool}"
    } else {
        println "Aucun outil de binning spécifié, utilisation de l'outil par défaut: metabat2"
    }

    // Définition des channels de base
    bam_merger_ch = ont_bam_ch.join(illumina_bam_ch)
    //bam_merger(bam_merger_ch)
    merged_bam_out = bam_merger(bam_merger_ch)
    //Channel metabat2_out_ch
    //Channel semibin2_out_ch
    //Channel comebin_out_ch
    // Logique de sélection de l'outil de binning
    switch (params.bintool) {
        // case 'metabat2':
        //     if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
        //         metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        //         metabat2_extra(metabat2_ch, extra_bam)
        //         metabat2_out = metabat2_extra.out
        //     }
        //     else {
        //         metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
        //         metabat2(metabat2_ch)
        //         metabat2_out = metabat2.out
        //     }
        //     break
        case 'metabat2':
            if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                metabat2_tmp_ch = merged_bam_out.join(extra_bam)
                bam_merger(metabat2_tmp_ch)
                metabat2_extra(assembly_ch, bam_merger.out)
                metabat2_out_ch = metabat2_extra.out
            }
            else {
                //metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                metabat2(assembly_ch.join(merged_bam_out))
                metabat2_out_ch = metabat2.out
            }
            break

        case 'semibin2':
        //possiblité de merger les extra reads
            semibin2_ch = assembly_ch.join(merged_bam_out)
            semibin2(semibin2_ch)
            semibin2_out_ch = semibin2.out
            break

        case 'comebin':
            comebin_ch = assembly_ch.join(merged_bam_out)
            comebin(comebin_ch)
            comebin_out_ch = comebin.out
            break

        default:
            println "L'outil spécifié (${params.bintool}) n'est pas reconnu. Utilisation de l'outil par défaut: metabat2"
            if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                metabat2_tmp_ch = merged_bam_out.join(extra_bam)
                bam_merger(metabat2_tmp_ch)
                metabat2_extra(assembly_ch, bam_merger.out)
                metabat2_out_ch = metabat2_extra.out
            }
            else {
                //metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                metabat2(assembly_ch, merged_bam_out)
                metabat2_out_ch = metabat2.out
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