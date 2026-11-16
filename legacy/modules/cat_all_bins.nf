process cat_all_bins {
    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(bins)
    output:
    tuple val(name) , path("all_bins.fa") 
    script:
    """
    cat ${bins}/bin.*.fa > all_bins.fa
    """
}