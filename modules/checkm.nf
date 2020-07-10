process checkm {
    maxForks 1
    label 'checkm'
    publishDir "${params.output}/${name}/classify/checkm/", mode: 'copy', pattern: "summary.txt"
    publishDir "${params.output}/${name}/classify/checkm/", mode: 'copy', pattern: "taxonomy.txt"
    publishDir "${params.output}/${name}/classify/checkm/", mode: 'copy', pattern: "*_checkm"
    publishDir "${params.output}/${name}/classify/checkm/", mode: 'copy', pattern: "*_checkm_plot"
    errorStrategy { task.exitStatus in 14..14 ? 'retry' : 'finish'}
    maxRetries 3 
    input:
    tuple val(name), path(bins_assemblies)
    output:
    tuple val(name), path("summary.txt")
    tuple path("${name}_checkm"), path("${name}_checkm_plot"), path("taxonomy.txt")
    
    script:
    """
    mkdir temporary
    mkdir ${name}_bin
    mv *.fa ${name}_bin/
    checkm lineage_wf --tmpdir temporary --pplacer_threads 4 -t ${task.cpus} --reduced_tree -x fa ${name}_bin ${name}_checkm > summary.txt
    checkm bin_qa_plot --image_type png -x fa ${name}_checkm ${name}_bin ${name}_checkm_plot
    checkm tree_qa ${name}_checkm > taxonomy.txt
     """
}

// checkm module is not use in the script at the moment but it is used in metawrap
// this module can be added for an additional check by the user just call it in the main script and input a channel outputted from a binning step