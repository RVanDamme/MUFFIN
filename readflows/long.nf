workflow [long_read_workflow]{
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

    Channel ont_input_ch
    // Assemblage avec des reads longs (ONT)
    if (!params.ont) error "ONT reads path must be specified for 'long' read type."
    ont_input_ch = Channel.fromPath("${params.ont}/*.fastq{,.gz}", checkIfExists: true).map { file -> tuple(file.baseName, file) }
    if (!params.skip_ont_qc) {
        ont_input_ch = ont_input_ch.flatMap { chopper(it) }
    }

    Channel assembly_ch

    // Assemblage avec des reads longs (ONT) utilisant Flye
    flye(ont_input_ch)
    assembly_ch = flye.out

    assembly_ch.flatMap { contigs ->
        minimap_polish(contigs, ont_input_ch)
        }.flatMap { polished ->
            racon(polished)
        }.flatMap { racon_out ->
            medaka(racon_out)
        }.flatMap { medaka_out ->
            pilong(medaka_out, ont_input_ch, params.polish_iteration)
        }.set { assembly_ch }
        //ajouter pilon

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

    //***************************************************
    // Assembly quality control
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
            Channel semibin2_ch = assembly_ch.join(ont_bam_ch)
            semibin2(semibin2_ch)
            Channel semibin2_out = semibin2.out
            break

        case 'comebin':
            Channel comebin_ch = assembly_ch.join(ont_bam_ch)
            comebin(comebin_ch)
            Channel comebin_out = comebin.out
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
}



