workflow [hybrid_workflow]{
    // Initialisation des variables pour les chemins des bases de données
    Channel<String> database_sourmash
    Channel<String> checkm_db_path

    // Configuration de la base de données Sourmash
    database_sourmash = params.sourmash_db ?
                        Channel.value(file(params.sourmash_db).toString()) :
                        sourmash_download_db().map { it.toString() }

    // Configuration de la base de données CheckM
    if (workflow.profile.contains('conda')) {
        // Utilisation de Conda, nécessite une configuration spécifique
        if (params.checkm_db) {
            // Chemin vers la base de données CheckM fourni
            checkm_db_path = checkm_setup_db(params.checkm_db, true).map { it.toString() }
        } else if (params.checkm_tar_db) {
            // Chemin vers l'archive tar de la base de données CheckM fourni
            checkm_db_path = checkm_setup_db(params.checkm_tar_db, false).map { it.toString() }
        } else {
            // Téléchargement et configuration de la base de données CheckM
            checkm_db_path = checkm_download_db()
                             .flatMap { checkm_setup_db(it, false) }
                             .map { it.toString() }
        }
    } else {
        // Utilisation d'un autre profil (ex. Docker), configuration par défaut
        checkm_db_path = Channel.value("/checkm_database/path.txt")
    }


    Channel illumina_input_ch, ont_input_ch

    if (!params.ont || !params.illumina) error "Both ONT and Illumina reads paths must be specified for 'hybride' read type."
    ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
    illumina_input_ch = Channel.fromFilePairs("${params.illumina}/*_R{1,2}.fastq{,.gz}", checkIfExists: true)
    if (!params.skip_ont_qc) {
        ont_input_ch = ont_input_ch.flatMap { chopper(it) }
    }
    if (!params.skip_ill_qc) {
        illumina_input_ch = illumina_input_ch.flatMap { fastp(it) }
    }


    Channel assembly_ch

    switch (params.assembler) {
        case "metaspades":
            Channel spades_ch = illumina_input_ch.join(ont_input_ch)

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
    Channel ont_bam_ch = assembly_ch
        .join(ont_input_ch)
        .map { assembly, reads -> [assembly, reads] }
        .flatMap { minimap2(it) }

    // Mapping additional ONT reads if specified
    if (params.extra_ont) {
        Channel ont_extra_bam_ch = assembly_ch
            .join(Channel.fromPath(params.extra_ont))
            .map { assembly, extraReads -> [assembly, extraReads] }
            .flatMap { extra_minimap2(it) }
    }

    // Mapping with BWA for Illumina reads
    Channel illumina_bam_ch = assembly_ch
        .join(illumina_input_ch)
        .map { assembly, reads -> [assembly, reads] }
        .flatMap { bwa(it) }

    // Mapping additional Illumina reads if specified
    if (params.extra_ill) {
        Channel illumina_extra_bam_ch = assembly_ch
            .join(Channel.fromPath(params.extra_ill))
            .map { assembly, extraReads -> [assembly, extraReads] }
            .flatMap { extra_bwa(it) }
    }


    //***************************************************
    // Assembly quality control-> may be pushed in a function
    //***************************************************
    if (params.reference) {
        //Channel ref_ch = params.reference
        Channel ref_ch = Channel.fromPath(params.reference)
        metaquast(assembly_ch, ref_ch)
        Channel metaquast_out_ch = metaquast.out
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
    Channel bam_merger_ch = ont_bam_ch.join(illumina_bam_ch)
    bam_merger(bam_merger_ch)
    Channel merged_bam_out = bam_merger.out

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
                metabat2_out = metabat2_extra.out
            }
            else {
                //metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                metabat2(assembly_ch, merged_bam_out)
                metabat2_out = metabat2.out
            }
            break

        case 'semibin2':
            Channel semibin2_ch = assembly_ch.join(merged_bam_out)
            semibin2(semibin2_ch)
            semibin2_out = semibin2.out
            break

        case 'comebin':
            Channel comebin_ch = assembly_ch.join(merged_bam_out)
            comebin(comebin_ch)
            comebin_out = comebin.out
            break

        default:
            println "L'outil spécifié (${params.bintool}) n'est pas reconnu. Utilisation de l'outil par défaut: metabat2"
            if (params.extra_ont || params.extra_ill ) { // check if differential coverage binning possible
                metabat2_tmp_ch = merged_bam_out.join(extra_bam)
                bam_merger(metabat2_tmp_ch)
                metabat2_extra(assembly_ch, bam_merger.out)
                metabat2_out = metabat2_extra.out
            }
            else {
                //metabat2_ch = assembly_ch.join(ont_bam_ch).join(illumina_bam_ch)
                metabat2(assembly_ch, merged_bam_out)
                metabat2_out = metabat2.out
            }
    }
}