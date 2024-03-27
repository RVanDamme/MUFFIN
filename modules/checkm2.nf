process checkm2 {
    maxForks 1
    label 'checkm2'

    conda 'bioconda::checkm2=1.0.1'
    
    publishDir "${params.output}/${name}/classify/checkm2/", mode: 'copy', pattern: "checkm2_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(bins_assemblies)
    val(bool)
    output:
    tuple val(name), path("checkm2_dir")
    
    script:
    // Détermine si params.checkm2db est fourni et ajoute --database_path à la commande si nécessaire
    def dbPathCmd = params.checkm2db ? "--database_path ${params.checkm2db}" : ""
    def lowmemCmd = params.checkm2_low ? "--lowmem" : ""
    
    """
    checkm2 predict --threads ${task.cpus} --input ${bins_assemblies}/metabat_bins.10.fa ${dbPathCmd} ${lowmemCmd} --output-directory checkm2_dir/ --database_path /home/nestpasici
    """
}

// Low memory mode
// If you are running CheckM2 on a device with limited RAM, you can use the 
// --lowmem option to reduce DIAMOND RAM use by half at the expense of longer runtime.
//--database_path /proj/cloacimonetes/NOBACKUP/Arnaud2024/spades_way/uniref100.KO.1.dmnd