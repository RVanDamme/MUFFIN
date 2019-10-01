process checkm {
    label 'checkm'
    publishDir "${params.output}/${name}_checkm_bins/", mode: 'copy', pattern: "*"
    input:
    set val(name), file(bins)
    output:
    //TBD
    script:
    """
    checkm lineage_wf -t ${task.cpus} -x fa ${bins} ${bins}"_checkm"
    checkm qa ${bins}"_checkm/lineage.ms" ${bins}"_checkm"
    checkm bin_qa_plot -x fa ${bins}"_checkm" ${bins} ${bins}"_checkm_plot"
    """
}