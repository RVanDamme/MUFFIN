process contig_list {
    label 'ubuntu'
    input:
    set val(name), file(bins)
    output:
    set val(name), file("ids/*.list")
    shell:\
    """
    mkdir ids
    for bin in \$(echo !{bins} |sed 's/[][]//g' | tr "," "\n")
        do
        bin_name=\$(basename \$bin )
        cat \$bin | grep -o -E "^>\\w+\\.\\w+" |sed 's/>//g'| tr -d "@" > ids/\$bin_name.contigs.list ;
        done ;
    """
}

