process refine2 {
    label 'ubuntu'
    input:
    output:
    tuple val("subset"), path("subset/substet_ill/subset_R{1,2}.fastq")
    tuple val("subset"), path("subset/subset_ont/subset.fastq")
    tuple val("subset"), path("subset/subset_rna/subset_R{1,2}.fastq.gz")
    shell:
    """
    wget https://osf.io/9xmh4/download -O subset_data.tar.gz
    tar -xzvf subset_data.tar.gz
    """
}