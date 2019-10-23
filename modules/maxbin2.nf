process maxbin2 {
    label 'maxbin2'
    if (params.out_maxbin == true ) { publishDir "${params.output}/${name}_maxbin2/", mode: 'copy', pattern: "maxbin_bin/*" }
    input:
    set val(name), file(assembly), file(ont), file(illumina)
    output:
    set val(name), file("maxbin_bin/")
    script:
    """
    run_MaxBin.pl -contig ${assembly}  -reads ${illumina[0]} -reads2 ${illumina[1]} -reads3 ${ont}  -out maxbin2 -thread ${task.cpus}
    mkdir maxbin_bin
    mv maxbin2.*.fa maxbin_bin/
    """
}  // add -prob_threshold 0.5 -markerset 40 ??