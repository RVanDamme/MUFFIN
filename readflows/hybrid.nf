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
    include {pilong} from '../modules/polish' params(output : params.output)
    include {minimap2} from '../modules/minimap2' //mapping for the binning 
    include {extra_minimap2} from '../modules/minimap2'
    include {bwa} from '../modules/bwa' //mapping for the binning
    include {bwa_bin} from '../modules/bwa' //mapping for the binning
    include {extra_bwa} from '../modules/bwa'
    //include {metabat2_extra} from '../modules/metabat2' params(output : params.output)    
    include {metabat2} from '../modules/metabat2' params(output : params.output)
    include {semibin2} from '../modules/semibin2' params(output : params.output, bining_model : params.bining_model, environment : params.environment)
    include {comebin} from '../modules/comebin' params(output : params.output)
    include {separateBins} from '../modules/bins_tools' params(output : params.output)
    include {bin_merger} from '../modules/bins_tools' params(output : params.output)
    include {bam_merger} from '../modules/samtools_merger' params(output : params.output)
    //include {contig_list} from '../modules/list_ids'
    include {cat_all_bins} from '../modules/cat_all_bins'
    //include {bwa_bin} from '../modules/bwa'  
    include {minimap2_bin} from '../modules/minimap2'
    include {metaquast} from '../modules/quast' params(output : params.output)
    include {unmapped_illumina_retrieve} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    include {unmapped_ont_retrieve} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    // include {illumina_reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    // include {ont_reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){
    include {sourmash_download_db} from '../modules/sourmashgetdatabase'
    include {checkm_download_db} from '../modules/checkmgetdatabases'
}
if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot") {
    //include {reads_retrieval} from '../modules/seqtk_retrieve_reads' params(output : params.output)
    include {checkm2} from '../modules/checkm2' params(output: params.output, checkm2_low: params.checkm2_low, checkm2db: params.checkm2db)
    include {sourmash_bins} from '../modules/sourmash' params(output: params.output, bintool : params.bintool)
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

include {readme_output} from '../modules/readme_output' params(output: params.output)
include {test} from '../modules/test_data_dll'

params.db_path = '/home/user/databases'
params.db_file = 'uniref100.KO.1.dmnd'

workflow hybrid_workflow{
    // Initialisation des variables pour les chemins des bases de données
    // if (params.modular=="full" | params.modular=="classify" | params.modular=="assem-class" | params.modular=="class-annot"){

    //     // Configuration de la base de données Sourmash
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
    

    //Channel illumina_input_ch, ont_input_ch

    if (!params.ont || !params.illumina) error "Both ONT and Illumina reads paths must be specified for 'hybride' read type."

    //ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
    ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map {file -> tuple(file.simpleName, file) }.view()

    illumina_input_ch = Channel.fromFilePairs("${params.illumina}/*_R{1,2}.fastq{,.gz}", checkIfExists: true)

    if (!params.skip_ont_qc) {
        ont_input_ch = chopper(ont_input_ch)
        //ont_input_ch = ont_input_ch.flatMap { chopper(it) }
    }
    if (!params.skip_ill_qc) {
        //illumina_input_ch = illumina_input_ch.flatMap { fastp(it) }
        illumina_input_ch = fastp(illumina_input_ch)
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
            //assembly_ch = flye.out
            // Chaîne de polissage simplifiée
            minimap_polish_ch = minimap_polish(flye.out.join(ont_input_ch))
            racon_ch = racon(ont_input_ch.join(flye.out).join(minimap_polish_ch))
            medaka_ch = medaka(racon_ch)
            assembly_ch = medaka_ch
            //assembly_ch = pilon(medaka_ch.join(illumina_input_ch), params.polish_iteration)

            // assembly_ch.flatMap { contigs ->
            //     minimap_polish(contigs, ont_input_ch)
            // }.flatMap { polished ->
            //     racon(polished)
            // }.flatMap { racon_out ->
            //     medaka(racon_out)
            // }.flatMap { medaka_out ->
            //     pilon(medaka_out, illumina_input_ch, params.polish_iteration)
            // }.set { assembly_ch }
            break

        default:
            error "Unrecognized assembler: ${params.assembler}. Should be 'metaspades' or 'metaflye'."
    }


    //*********
    // Mapping
    //*********
    
    // Mapping with Minimap2 for ONT reads
    ont_bam_ch = minimap2(assembly_ch.join(ont_input_ch))

    // Mapping additional ONT reads if specified
    if (params.extra_ont) {
        ont_extra_bam_ch = minimap2(assembly_ch.join(extra_ont_ch))
    }

    // Mapping with BWA for Illumina reads
    illumina_bam_ch = bwa(assembly_ch.join(illumina_input_ch))

    // Mapping additional Illumina reads if specified
    if (params.extra_ill) {
        illumina_extra_bam_ch = bwa(assembly_ch.join(params.extra_ill))
    }


    //***************************************************
    // Assembly quality control-> may be pushed in a function
    //***************************************************
    if (params.reference) {
        //Channel ref_ch = params.reference
        ref_ch = Channel.fromPath(params.reference)
        // metaquast(assembly_ch, ref_ch)
        // metaquast_out_ch = metaquast.out
        metaquast_out_ch = metaquast(assembly_ch.join(ref_ch))
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
    merged_bam_out = bam_merger(bam_merger_ch)
    // Logique de sélection de l'outil de binning
    switch (params.bintool) {
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
        //else {classify_ch = metabat2_out_ch.join(semibin2_out_ch).join(comebin_out_ch)}
        else{
            switch (params.bintool) {
                case 'metabat2':
                    bin_ch = metabat2.out
                    break

                case 'semibin2':
                    bin_ch = semibin2.out
                    break

                case 'comebin':
                    bin_ch = comebin.out
                    break

                default:
                    bin_ch = metabat2.out

            }
            classify_ch = bin_ch
        }


        // Configuration de la base de données Sourmash
        if (params.sourmash_db) { database_sourmash = file(params.sourmash_db) }
        else {
            sourmash_download_db() 
            database_sourmash = sourmash_download_db.out
        }   
        if (!params.checkm2db){
            // Vérification de l'existence du dossier et du fichier
            path_exists = file(params.db_path).exists() && file("${params.db_path}/${params.db_file}").exists()
            //println "checkm2 db install process"
            //println path_exists
            // Si le chemin n'existe pas ou si le fichier n'est pas trouvé.
            if( !path_exists | params.checkm2db_force_update) {
                checkm_download_db()
            } else {
                println "Le dossier et le fichier spécifié existent déjà."
            }
        }

    
        //*************************
        // Bins classify workflow
        //*************************

        //checkm of the final assemblies
        //checkm(classify_ch.groupTuple(by:0)) //checkm QC of the bins
        checkm2(classify_ch, checkm_download_db.out)
        // classify_ch.view()

        if (!params.skip_bin_sorting){
            separateBins(checkm2.out.join(classify_ch))
            classify_ch = separateBins.out[0]
            bad_bins_ch = separateBins.out[1]
            //bad_bins_ch.view()
            
            if (!params.skip_bad_reads_recovery){
                merged_bin_ch = bin_merger(classify_ch)

                bwa_bin_ch = bwa_bin(merged_bin_ch.join(illumina_input_ch))
                minimap_bin_ch = minimap2_bin(merged_bin_ch.join(ont_input_ch))

                // ont_input_ch = ont_reads_retrieval(minimap_bin_ch.join(ont_input_ch))
                // illumina_input_ch = illumina_reads_retrieval(bwa_bin_ch.join(illumina_input_ch))

                unmapped_illumina_retrieve(bwa_bin_ch.join(illumina_input_ch))
                unmapped_ont_retrieve(minimap_bin_ch.join(ont_input_ch))
                //bwa_bin(bins_ready_ch, illumina_input_ch)
            }
        }
        
        //checkm2_out_ch = checkm2.out 
        classify_ch.flatMap { name, paths ->
            paths.collect { path -> tuple(name, path) }
        }
        .set { bins_ready_ch }

        if (!params.skip_pilon && params.assembler == 'metaflye'){
            pilon(bins_ready_ch.join(illumina_input_ch), params.polish_iteration)
        }

        


        //sourmash classification using gtdb database
        //sourmash_bins(classify_ch,database_sourmash) // fast classification using sourmash with the gtdb (not the best classification but really fast and good for primarly result)
        //sourmash_checkm_parser(checkm.out[0],sourmash_bins.out.collect()) //parsing the result of sourmash and checkm in a single result file

        //sourmash_bins(bins_ready_ch,database_sourmash)
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
            switch (params.bintool) {
                case 'metabat2':
                    bin_ch = metabat2.out
                    break

                case 'semibin2':
                    bin_ch = semibin2.out
                    break

                case 'comebin':
                    bin_ch = comebin.out
                    break

                default:
                    bin_ch = metabat2.out

            }
            bins_input_ch = bin_ch 
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
        // bins_input_ch.flatMap { name, paths ->
        //     paths.collect { path -> tuple(name, path) }
        // }
        // .set { bins_input_ready_ch }

        // eggnog_bin_ch = bins_input_ready_ch.combine(eggnog_db)
        // eggnog_bin(eggnog_bin_ch) //annotate the bins
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