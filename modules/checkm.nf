process checkm {
    maxForks 1
    label 'checkm'
    publishDir "${params.output}/${name}_checkm_bins/", mode: 'copy', pattern: "*"
    input:
    set val(name), val(bin_id), file(bins_assemblies)
    output:
    //TBD
    script:
    """
    mkdir ${bin_id}_bin
    mv *.fa ${bin_id}_bin/
    checkm lineage_wf -t ${task.cpus} -x fa ${bin_id}_bin ${bin_id}_checkm > summary.txt
    checkm qa ${bin_id}_checkm/lineage.ms ${bin_id}_checkm
    checkm bin_qa_plot -x fa ${bin_id}_checkm ${bin_id}_bin ${bin_id}_checkm_plot
    """
}
// checkm module is not use in the script at the moment but it is used in metawrap
// this module can be added for an additional check by the user just call it in the main script and input a channel outputted from a binning step