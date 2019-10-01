process metabat2 {
    label 'metabat2'
    input:
    set val(name), file(assembly), file(ont_bam), file(illumina_bam)
    output:
    set val(name), file("bins_dir/metabat_bin")
    script:
    """
    metabat -i ${assembly} ${ont_bam} ${illumina_bam} -o bins_dir/metabat_bin -t ${task.cpus}
    """
}