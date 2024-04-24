// process comebin {
//     maxForks 1
//     label 'comebin'

//     conda 'bioconda::COMEBin=1.0.3'

//     publishDir "${params.output}/${name}/assemble/binning/comebin/", mode: 'copy', pattern: "bins_dir"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     //tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
//     tuple val(name), path(assembly), path(bam_files)
//     output:
//     tuple val(name), path("bins_dir")
    
//     script:
//     """
//     #!/bin/bash
//     # Etape 1: Extraire les longueurs des contigs et les trier
//     readarray -t lengths < <(awk -F'_' '/^>/ {print \$4}' "${assembly}" | sort -nr)

//     # Etape 2: Calculer la somme totale des longueurs
//     total_length=0
//     for length in "\${lengths[@]}"; do
//         ((total_length += length))
//     done

//     # Etape 3: Trouver le N50
//     half_total=\$((total_length / 2))
//     cumulative_length=0
//     N50=0
//     for length in "\${lengths[@]}"; do
//         ((cumulative_length += length))
//         if [[ cumulative_length -ge half_total ]]; then
//             N50=\$length
//             break
//         fi
//     done

//     # Definir la temperature dans la fonction de perte en fonction du N50
//     loss_temp=\$(awk -v n50=\$N50 'BEGIN{print (n50 > 10000) ? 0.07 : 0.15}')
//     echo "loss temp"
//     echo \$loss_temp
    
//     temp_bam_dir="./temp_bam_dir"
//     mkdir -p \$temp_bam_dir

//     # Copie du fichier BAM dans le répertoire temporaire
//     cp ${bam_files} \$temp_bam_dir

//     echo ${task.cpus}
//     run_comebin.sh -t ${task.cpus} -a ${assembly} -o bins_dir/ -l \$loss_temp -p \$temp_bam_dir/
//     rm -r \$temp_bam_dir
//     """
// }


process comebin {
    maxForks 1
    label 'comebin'

    conda 'bioconda::COMEBin=1.0.3 bioconda::n50'

    publishDir "${params.output}/${name}/assemble/binning/comebin/", mode: 'copy', pattern: "bins_dir/comebin_res/comebin_res_bins"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    //tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
    tuple val(name), path(assembly), path(bam_files)
    output:
    tuple val(name), path("bins_dir/comebin_res/comebin_res_bins")
    
    script:
    """
    #!/bin/bash
    N50=\$(n50 ${assembly})
    echo "N50 calculé: \$N50"

    # Definir la temperature dans la fonction de perte en fonction du N50
    loss_temp=\$(awk -v n50=\$N50 'BEGIN{print (n50 > 10000) ? 0.07 : 0.15}')
    echo "loss temp"
    echo \$loss_temp
    
    temp_bam_dir="./temp_bam_dir"
    mkdir -p \$temp_bam_dir

    # Copie du fichier BAM dans le répertoire temporaire
    cp ${bam_files} \$temp_bam_dir

    echo ${task.cpus}
    run_comebin.sh -t ${task.cpus} -a ${assembly} -o bins_dir/ -l \$loss_temp -p \$temp_bam_dir/
    rm -r \$temp_bam_dir
    """
}

//run_comebin.sh -t ${task.cpus} -a ${assembly} -o bins_dir/comebin_bins -l \$loss_temp -p *.bam