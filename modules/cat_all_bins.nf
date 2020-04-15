process cat_all_bins {
    label 'ubuntu'
    input:
    tuple val(name), file(bins)
    output:
    tuple val(name) , file("all_bins.fa") 
    script:
    """
    cat ${bins}/bin.*.fa > all_bins.fa
    """
}