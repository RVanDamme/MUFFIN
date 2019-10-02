process refine2 {
    label 'metawrap'
    if (params.out_metawrap == true ) { publishDir "${params.output}/${name}_refined_bins/", mode: 'copy', pattern: "refined_bins/*" }
    input:
    set val(name1), file(bins1), file(bins2)
    output:
    set val(name1), file("refined_bins")
    script:
    """
    metaWRAP bin_refinement -t ${tasks.cpus} -m ${tasks.memory} -o refined_bins -A ${bins1} -B ${bins2} -o refined_bins
    """
}

process refine3 {
    label 'metawrap'
    if (params.out_metawrap == true ) { publishDir "${params.output}/${name}_refined_bins/", mode: 'copy', pattern: "refined_bins/*" }
    input:
    set val(name1), file(bins1), file(bins2), file(bins3)
    output:
    set val(name1), file("refined_bins")
    script:
    """
    metaWRAP bin_refinement -t ${tasks.cpus} -m ${tasks.memory} -o refined_bins -A ${bins1} -B ${bins2} -C ${bins3} -o refined_bins
    """
}