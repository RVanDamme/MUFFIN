process separateBins {

    label 'ubuntu'
    //publishDir "${params.output}/${name}/annotate/", mode: 'copy', pattern: "parser_result/*"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5

    input:
    tuple val(name), path(checkm2_res_file)
    tuple val(name), path(bins)

    output:
    tuple val(name), path(good_bins)
    tuple val(name), path(bad_bins)

    script:
    """
    # Initialiser des fichiers temporaires pour stocker les chemins
    good_bin_paths=\$(mktemp)
    bad_bin_paths=\$(mktemp)

    # Lire le fichier de données et séparer les chemins des bins
    awk -v good=\$good_bin_paths -v bad=\$bad_bin_paths -v bin_dir="${bins}" 'NR > 1 {
        if ((\$2 - 5*\$3) > 50)
            print bin_dir \$1 ".fa" >> good;
        else
            print bin_dir \$1 ".fa" >> bad;
    }' ${checkm2_res_file}

    # Utiliser 'cat' pour créer des listes de fichiers Nextflow à partir des chemins
    cat \$good_bin_paths > good_bins
    cat \$bad_bin_paths > bad_bins

    # Nettoyage
    rm -f \$good_bin_paths \$bad_bin_paths
    """
}