process refine2 {
    //label 'metawrap'
    conda '/home/renaud/miniconda3/envs/metawrap-env'
    if (params.out_metawrap == true ) { publishDir "${params.output}/${name}_refined_bins/", mode: 'copy', pattern: "refined_bins/*" }
    input:
    set val(name1), file(bins1), file(bins2)
    file(path)
    output:
    set val(name1), file("refined_bins/metawrap_70_10_bins/*.fa")
    shell:
    """
    mem=\$(echo !{task.memory} | sed 's/ GB//g')
    path_db=\$(cat !{path})
    echo \$path_db
    echo -e "\$path_db" | checkm data setRoot
    echo "checkm done"
    metawrap bin_refinement -t !{task.cpus} -m \$mem -o refined_bins -A !{bins1} -B !{bins2} -o refined_bins
    """
}

process refine3 {
    //label 'metawrap'
    label 'ubuntu'
    //label 'metawrap'
    conda '/home/renaud/miniconda3/envs/metawrap-env'
    if (params.out_metawrap == true ) { publishDir "${params.output}/${name}_refined_bins/", mode: 'copy', pattern: "refined_bins/*" }
    input:
    set val(name1), file(bins1), file(bins2), file(bins3)
    file(path)
    output:
    set val(name1), file("refined_bins/metawrap_70_10_bins/*.fa")
    shell:
    """
    mem=\$(echo !{task.memory} | sed 's/ GB//g')
    path_db=\$(cat !{path})
    echo \$path_db
    echo -e "\$path_db" | checkm data setRoot
    echo "checkm done"
    metawrap bin_refinement -t !{task.cpus} -m \$mem -o refined_bins -A !{bins1} -B !{bins2} -C !{bins3}

    """
}

    // path_db=\$(cat !{path})
    // echo \$path_db
    // echo -e "\$path_db" | checkm data setRoot
    // echo "checkm done"
    // metawrap bin_refinement -t !{task.cpus} -m !{task.memory} -o refined_bins -A !{bins1} -B !{bins2} -C !{bins3} -o refined_bins
