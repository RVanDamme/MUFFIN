process checkm {
    maxForks 1
    label 'checkm'
    publishDir "${params.output}/${name}/checkm_bins/", mode: 'copy', pattern: "summary.txt"
    publishDir "${params.output}/${name}/checkm_bins/", mode: 'copy', pattern: "taxonomy.txt"
    publishDir "${params.output}/${name}/checkm_bins/", mode: 'copy', pattern: "*_checkm"
    publishDir "${params.output}/${name}/checkm_bins/", mode: 'copy', pattern: "*_checkm_plot"
    input:
    set val(name), file(bins_assemblies)
    output:
    set val(name), file("summary.txt")
    set file("${name}_checkm"), file("${name}_checkm_plot"), file("taxonomy.txt")
    
    script:
    """
    mkdir temporary
    mkdir ${name}_bin
    mv *.fa ${name}_bin/
    checkm lineage_wf --tmpdir temporary --pplacer_threads ${task.cpus} -t ${task.cpus} --reduced_tree -x fa ${name}_bin ${name}_checkm > summary.txt
    checkm bin_qa_plot --image_type png -x fa ${name}_checkm ${name}_bin ${name}_checkm_plot
    checkm tree_qa ${name}_checkm > taxonomy.txt
     """
}

// checkm module is not use in the script at the moment but it is used in metawrap
// this module can be added for an additional check by the user just call it in the main script and input a channel outputted from a binning step