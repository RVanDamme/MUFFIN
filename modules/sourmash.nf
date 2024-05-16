process sourmash_genome_size { //deprecated since flye 2.8 update
    label 'sourmash' 

    conda 'bioconda::sourmash=2.0.1 '

    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(ont)
    path(json)
    output:
    tuple val(name), path(ont), path('genome_size.txt')
    shell:
    """
    echo "100M" >genome_size.txt
    """
/*
    sourmash compute -p !{task.cpus} --scaled 10000 -k 31 !{ont} -o !{name}.sig 
    sourmash lca gather  !{name}.sig !{json} --ignore-abundance -o metagenomic-composition.txt
    sum_ont=\$(cat metagenomic-composition.txt | cut -d ',' -f 1 | paste -sd+ | bc)
    total_m_ont=\$(bc -l <<< "scale=2 ; \$sum_ont /10^6")
    if (( \$(echo "\$total_m_ont < 100" |bc -l) ));
        then echo "100M" >genome_size.txt;
        else echo \$total_m_ont"M" >genome_size.txt;
    fi*/

}

process sourmash_bins {
    label 'sourmash' 

    conda 'bioconda::sourmash=2.0.1 '

    publishDir "${params.output}/${name}/classify/sourmash/", mode: 'copy', pattern: "*.txt"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(bins)
    path(json)
    output:
    path('*.txt')
    shell:
    """
    bin_id=\$(basename ${bins} | sed -r "s/\\.\\w+//2")
    sourmash compute -p ${task.cpus} --scaled 10000 -k 31 ${bins} -o ${bins}.sig
    sourmash lca classify --query ${bins}.sig --db ${json} > \$bin_id.txt   
    """
}

process sourmash_ill {
    label 'sourmash' 

    conda 'bioconda::sourmash=2.0.1 '

    publishDir "${params.output}/${name}/classify/not_mapped_sourmash/ill/", mode: 'copy', pattern: "*.kreport"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(illumina_reads)
    path(lineages)
    path(full_db)
    output:
    path('*.kreport')
    shell:
    """
    ill_forward=\$(basename ${illumina_reads[0]} | sed -r "s/\\.\\w+//2")
    ill_reverse=\$(basename ${illumina_reads[1]} | sed -r "s/\\.\\w+//2") 

    sourmash sketch dna -p k=31,scaled=10000 -o \$ill_forward.sig.zip ${illumina_reads[0]}
    sourmash gather \$ill_forward.sig.zip ${full_db} -k 31 -o \$ill_forward.gather.k31.csv
    sourmash tax metagenome --gather-csv \$ill_forward.gather.k31.csv --taxonomy ${lineages} --output-format kreport > \$ill_forward.kreport

    sourmash sketch dna -p k=31,scaled=10000 -o \$ill_reverse.sig.zip ${illumina_reads[1]}
    sourmash gather \$ill_reverse.sig.zip ${full_db} -k 31 -o \$ill_reverse.gather.k31.csv
    sourmash tax metagenome --gather-csv \$ill_reverse.gather.k31.csv --taxonomy ${lineages} --output-format kreport > \$ill_reverse.kreport

    """
}

process sourmash_ont {
    label 'sourmash' 

    conda 'bioconda::sourmash=2.0.1 '

    publishDir "${params.output}/${name}/classify/not_mapped_sourmash/ont/", mode: 'copy', pattern: "*.kreport"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(ont_reads)
    path(lineages)
    path(full_db)
    output:
    path('*.kreport')
    shell:
    """
    ont_id=\$(basename ${ont_reads} | sed -r "s/\\.\\w+//2")

    sourmash sketch dna -p k=31,scaled=5000 -o \$ont_id.sig.zip ${ont_reads}
    sourmash gather \$ont_id.sig.zip ${full_db} -k 31 -o \$ont_id.gather.k31.csv
    sourmash tax metagenome --gather-csv \$ont_id.gather.k31.csv --taxonomy ${lineages} --output-format kreport > \$ont_id.kreport 
    
    """
}


// process sourmash_bins {
//     label 'sourmash' 

//     conda 'bioconda::sourmash=2.0.1 '

//     publishDir "${params.output}/${name}/classify/sourmash/", mode: 'copy', pattern: "*.txt"
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     input:
//     tuple val(name), path(bins)
//     path(json)
//     output:
//     path('*.txt')
//     shell:
//     """
//     sourmash compute -p ${task.cpus} --scaled 10000 -k 31 ${bins}/* -o sourmash.sig
//     sourmash lca classify --query sourmash.sig --db ${json} > souremash_output.txt   
//     """
// }
