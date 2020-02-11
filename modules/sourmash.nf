process sourmash_genome_size {
    label 'sourmash' 
    input:
    set val(name), file(ont)
    file(json)
    output:
    set val(name), file(ont), file('genome_size.txt')
    shell:
    """
    sourmash compute -p !{task.cpus} --scaled 10000 -k 31 !{ont} -o !{name}.sig 
    sourmash lca gather  !{name}.sig !{json} --ignore-abundance -o metagenomic-composition.txt
    sum_ont=\$(cat metagenomic-composition.txt | cut -d ',' -f 1 | paste -sd+ | bc)
    total_m_ont=\$(bc -l <<< "scale=2 ; \$sum_ont /10^6")
    if (( \$(echo "\$total_m_ont < 100" |bc -l) ));
        then echo "100M" >genome_size.txt;
        else echo \$total_m_ont"M" >genome_size.txt;
    fi
    """
}

process sourmash_bins {
    label 'sourmash' 
    publishDir "${params.output}/${name}/sourmash/${bin_id}/", mode: 'copy', pattern: "*.txt"
    input:
    set val(name), file(bins)
    file(json)
    output:
    file('*.txt')
    shell:
    """
    bin_id=\$(basename !{bins} | sed -r "s/\\.\\w+//2")
    sourmash compute -p !{task.cpus} --scaled 10000 -k 31 !{bins} -o !{bins}.sig
    sourmash lca classify --query !{bins}.sig --db !{json} > \$bin_id.txt   
    """
}
