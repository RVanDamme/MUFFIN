// process comebin {
//     maxForks 1
//     label 'semibin2'
//     publishDir "\${params.output}/\${name}/assemble/binning/semibin2/", mode: 'copy', pattern: "bins_dir"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     //tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
//     tuple val(name), path(assembly), path(bam_files)
//     output:
//     tuple val(name), path("bins_dir")
    
//     script:
//     """
//     mkdir temp;
//     cat \${assembly} | awk '!/^>/ { print "%s", \$0; n = "\n" } /^>/ { print n \$0; n = ""} END { printf "%s", n }' | sed '/^>/ d'| awk '{ print length(\$0) }' | sort -gr > temp/contig_lengths.txt;

//     #number of contigs
//     Y=$(cat temp/contig_lengths.txt | wc -l);

//     #sum of contig_lengths
//     X=$(paste -sd+ temp/contig_lengths.txt | bc);

//     # cumulative contig lengths
//     awk 'BEGIN {sum=0} {sum= sum+\$0; print sum}' temp/contig_lengths.txt > temp/contig_lengths_cum.txt

//     # get cumulative contig contributions (%) to the entire assembly
//     awk -v var=$X 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) \$i/=var;}1' temp/contig_lengths_cum.txt > temp/cum_perc.txt; 

//     # join results
//     paste temp/contig_lengths.txt temp/cum_perc.txt > temp/matrix.txt;

//     N50=$(awk '$2 >= 0.50' temp/matrix.txt |head -1| awk '{ print $1}');
//     rm -r temp;

//     run_comebin.sh [options] -t \${task.cpus} -a \${assembly} -o bins_dir/comebin_bins -p \$ *.bam
//     """
// }

// process n50 {
//     maxForks 1
//     label 'semibin2'
//     publishDir "\${params.output}/\${name}/assemble/binning/semibin2/", mode: 'copy', pattern: "bins_dir"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     //tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
//     tuple val(name), path(assembly), path(bam_files)
//     output:
//     tuple val(name), path("n50")
    
//     script:
//     """
//     #!/bin/bash
//     assembly="\${assembly}"

//     # Calculer les longueurs de contigs et trier en ordre décroissant
//     contig_lengths=\$(awk '!/^>/ { printf "%s", \$0; next } /^>/ { if(NR > 1) print n; n=0 } { n += length(\$0) } END { print n }' \$assembly | sort -nr)

//     # Calculer la somme totale des longueurs de contigs
//     total_length=\$(echo "\$contig_lengths" | awk '{sum += \$1} END {print sum}')

//     # Calculer N50
//     cumulative_length=0
//     for length in \$contig_lengths; do
//         cumulative_length=\$((\$cumulative_length + \$length))
//         if [ \$cumulative_length -ge \$(echo \$total_length / 2 | bc) ]; then
//             N50=\$length
//             break
//         fi
//     done

//     echo "N50: \$N50" 
//     """
// }


process comebin {
    maxForks 1
    label 'comebin'

    conda 'bioconda::COMEBin=1.0.3'

    publishDir "${params.output}/${name}/assemble/binning/semibin2/", mode: 'copy', pattern: "bins_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    //tuple val(name), path(assembly), path(ont_bam).optional(), path(illumina_bam).optional(), path(extra_bam).optional()
    tuple val(name), path(assembly), path(bam_files)
    output:
    tuple val(name), path("bins_dir")
    
    script:
    """
    #!/bin/bash
    assembly="${assembly}"

    # Calculer les longueurs de contigs et trier en ordre décroissant
    contig_lengths=\$(awk '!/^>/ { printf "%s", \$0; next } /^>/ { if(NR > 1) print n; n=0 } { n += length(\$0) } END { print n }' \$assembly | sort -nr)

    # Calculer la somme totale des longueurs de contigs
    total_length=\$(echo "\$contig_lengths" | awk '{sum += \$1} END {print sum}')

    # Calculer N50
    cumulative_length=0
    for length in \$contig_lengths; do
        cumulative_length=\$((\$cumulative_length + \$length))
        if [ \$cumulative_length -ge \$(echo \$total_length / 2 | bc) ]; then
            N50=\$length
            break
        fi
    done

    echo "N50: \$N50"

    # Définir la température dans la fonction de perte en fonction du N50
    loss_temp=\$(awk -v n50=\$N50 'BEGIN{print (n50 > 10000) ? 0.07 : 0.15}')

    run_comebin.sh -t ${task.cpus} -a ${assembly} -o bins_dir/comebin_bins -l \$loss_temp -p *.bam
    """
}