process refine2 {
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5   
    // uncomment the if else if you encounter issue using conda with metawrap
    //if (workflow.profile.contains('conda')) {conda '/path/to/miniconda3/envs/metawrap-env'}
    //else {label 'metawrap'}
    label 'metawrap'
    publishDir "${params.output}/${name}/assemble/binning/metawrap_refined_bins/", mode: 'copy', pattern: "metawrap_bins/*" 
    publishDir "${params.output}/${name}/assemble/binning/metawrap_refined_bins/", mode: 'copy', pattern: "${name}_binning_stats.txt" 
    input:
    tuple val(name), path(bins1), path(bins2)
    path(path)
    output:
    tuple val(name), path("metawrap_bins/*.fa")
    path("${name}_binning_stats.txt")
    shell:
    """
    mem=\$(echo ${task.memory} | sed 's/g//g')
    path_db=\$(cat ${path})
    echo \$path_db
    checkm data setRoot \$path_db
    echo "checkm done"
    metawrap bin_refinement -o refined_bins -A ${bins1} -B ${bins2} -o refined_bins -t ${task.cpus} -m \$mem 
    mkdir metawrap_bins/
    mv refined_bins/metawrap_70_10_bins/*.fa metawrap_bins/
    mv refined_bins/metawrap_70_10_bins.stats ${name}_binning_stats.txt
    """

}

process refine3 {
    // uncomment the if else if you encounter issue using conda with metawrap
    //if (workflow.profile.contains('conda')) {conda '/path/to/miniconda3/envs/metawrap-env'}
    //else {label 'metawrap'}
    label 'metawrap'
    publishDir "${params.output}/${name}/assemble/binning/metawrap_refined_bins/", mode: 'copy', pattern: "metawrap_bins/*" 
    publishDir "${params.output}/${name}/assemble/binning/metawrap_refined_bins/", mode: 'copy', pattern: "${name}_binning_stats.txt" 
    errorStrategy { task.exitStatus in 1..1 ? 'retry' : 'finish'}
    maxRetries 1
    input:
        tuple val(name), path(bins1), path(bins2), path(bins3)
        path(path)
    output:
    tuple val(name), path("metawrap_bins/*.fa")
    path("${name}_binning_stats.txt")
    shell:
    if (task.attempt == 1)
    """
    mem=\$(echo ${task.memory} | sed 's/g//g')
    path_db=\$(cat ${path})
    echo \$path_db
    checkm data setRoot \$path_db
    echo "checkm done"
    metawrap bin_refinement -o refined_bins -A ${bins2} -B ${bins3} -C ${bins1} -t ${task.cpus} -m \$mem 
    mkdir metawrap_bins/
    mv refined_bins/metawrap_70_10_bins/*.fa metawrap_bins/
    mv refined_bins/metawrap_70_10_bins.stats ${name}_binning_stats.txt
    """
    else if (task.attempt == 2)
    """
    mem=\$(echo ${task.memory} | sed 's/g//g')
    path_db=\$(cat ${path})
    echo \$path_db
    echo -e "cat << EOF\\n\$path\\nEOF\\n" | checkm data setRoot 
    echo "checkm done"
    metawrap bin_refinement -o refined_bins -A ${bins2} -B ${bins1} -t ${task.cpus} -m \$mem 
    mkdir metawrap_bins/
    mv refined_bins/metawrap_70_10_bins/*.fa metawrap_bins/
    mv refined_bins/metawrap_70_10_bins.stats ${name}_binning_stats.txt
    """
    else 
    error "please pick the bins you want and submit them in the classification step then select the bins for the annotation step"
}

    // cp -r /home/renaud/mafin_modul/metawrap/metawrap_bins .
    // cp /home/renaud/mafin_modul/metawrap/S_41_17_Cf_binning_stats.txt .

process norefine {     
    label 'ubuntu'
    input:
    tuple val(name), path(bins)
    output:
    tuple val(name), path("${bins}/*.fa")
    shell:\
    """
    mkdir norefine
    cp ${bins}/* norefine/
    """
}
