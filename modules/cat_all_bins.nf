process cat_all_bins {
    //label 'ubuntu'
    input:
    set val(name), file(bins)
    output:
    set val(name) , file("all_bins.fa") emit: cat_bins
    set val(name), file("${bin}/bin*/") emit: independent_bin
    script:
    """
    cat ${bins}/*/*.{fa,fasta} > all_bins.fa
    """
}