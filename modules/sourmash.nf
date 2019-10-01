process sourmash {
    label 'sourmash'
    input:
    set val(name), file(ont)
    file(database)
    output:
    set val(name), file(ont), val('genome_size.txt')
    shell:
    """
    sourmash compute -p !{task.cpus} --scaled 10000 -k 31 !{ont} -o !{name}.sig 
    sourmash lca gather  !{name}.sig !{database} --ignore-abundance -o metagenomic-composition.txt
    sum_!{name}=\$(cat metagenomic-composition.txt | cut -d ',' -f 1 | paste -sd+ | bc)
    total_m_!{name}=\$(bc -l <<< "scale=2 ; \$sum_!{name} /10^6")
    echo \$total_m_!{name}"M" >genome_size.txt
    """
}