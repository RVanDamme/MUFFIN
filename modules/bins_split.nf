process separateBins {

    label 'ubuntu'
    publishDir "${params.output}/${name}/classify/sorted_bins/", mode: 'copy', pattern: "*bin_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5

    input:
    tuple val(name), path(checkm2_res_file)
    tuple val(name), path(bins_dir)

    output:
    tuple val(name), path("good_bin_dir/*.fa")
    tuple val(name), path("bad_bin_dir/*.fa")

    script:
    """
    # Initialiser des fichiers temporaires pour stocker les chemins

    echo ${checkm2_res_file}
    echo ${bins_dir}

    good_bin_dir="./good_bin_dir/"
    bad_bin_dir="./bad_bin_dir/"

    mkdir -p "\$good_bin_dir"
    mkdir -p "\$bad_bin_dir"
    
    awk -v dir="${bins_dir}" -v good_dir="\$good_bin_dir" -v bad_dir="\$bad_bin_dir" 'NR > 1 {
    if ((\$2 - 5*\$3) > 50)
        print "cp " ${bins_dir} \$1 ".fa " good_dir;
    else
        print \$1;
    }' "${checkm2_res_file}"

    """
}

// awk -v dir="${bins_dir}" -v good_dir="\$good_bin_dir" -v bad_dir="\$bad_bin_dir" 'NR > 1 {
//     if ((\$2 - 5*\$3) > 50)
//         system("cp " dir \$1 ".fa " good_dir);
//     else
//         system("cp " dir \$1 ".fa " bad_dir);
//     }' "${checkm2_res_file}"