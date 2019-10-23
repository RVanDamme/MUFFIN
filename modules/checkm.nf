process checkm {
    label 'checkm'
    publishDir "${params.output}/${name}_checkm_bins/", mode: 'copy', pattern: "*"
    input:
    set val(name), file(bins_assemblies)
    output:
    //TBD
    script:
    """
    mkdir bins
    mv *_polished.fasta bins/
    checkm lineage_wf -t ${task.cpus} -x fasta bins bins_checkm
    checkm qa bins_checkm/lineage.ms bins_checkm
    checkm bin_qa_plot -x fasta bins_checkm bins bins_checkm_plot
    """
}
// checkm module is not use in the script at the moment but it is used in metawrap
// this module can be added for an additional check by the user just call it in the main script and input a channel outputted from a binning step