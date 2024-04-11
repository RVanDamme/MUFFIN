process separateBins {

    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5

    input:
    tuple val(name), path(checkm2_dir), path(bins_dir)

    output:
    tuple val(name), path("good_bin_dir/*.fa")
    tuple val(name), path("bad_bin_dir/*.fa")

    script:
    """
    good_bin_dir="./good_bin_dir/"
    bad_bin_dir="./bad_bin_dir/"

    mkdir -p "\$good_bin_dir"
    mkdir -p "\$bad_bin_dir"
    
    awk -v dir="${bins_dir}/" -v good_dir="\$good_bin_dir" -v bad_dir="\$bad_bin_dir" 'NR > 1 {
    if ((\$2 - 5*\$3) > 50)
        system("cp " dir \$1 ".fa " good_dir);
    else
        system("cp " dir \$1 ".fa " bad_dir);
    }' "${checkm2_dir}/quality_report.tsv"

    cp -r \$good_bin_dir "${params.output}/${name}/classify/sorted_bins/"
    cp -r \$bad_bin_dir "${params.output}/${name}/classify/sorted_bins/"

    """
}

process bin_filter {

    label 'ubuntu'
    publishDir "${params.output}/${name}/classify/sorted_bins/bin_filtered", mode: 'copy', pattern: "good_bin_dir/*.fa"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5

    input:
    tuple val(name), path(checkm2_dir), path(bins_dir)

    output:
    tuple val(name), path("good_bin_dir/*.fa")

    script:
    """
    good_bin_dir="./good_bin_dir/"

    mkdir -p "\$good_bin_dir"
    
    awk -v dir="${bins_dir}/" -v good_dir="\$good_bin_dir" 'NR > 1 {
    if ((\$2 - 5*\$3) > 50)
        system("cp " dir \$1 ".fa " good_dir);

    }' "${checkm2_dir}/quality_report.tsv"

    """
}

process bin_merger {

    label 'ubuntu'
    publishDir "${params.output}/${name}/classify/sorted_bins/", mode: 'copy', pattern: "merged_bin_assembly.fa"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate'}
    maxRetries = 5

    input:
    tuple val(name), path(bins)

    output:
    tuple val(name), path("merged_bin_assembly.fa")

    script:
    """
    cat ${bins} > merged_bin_assembly.fa
    """
}