process test {
    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    output:
    tuple val("subset"), path("ill/*")
    tuple val("subset"), path("ont/subset.fastq")
    tuple val("subset"), path("rna/*")
    shell:
    """
    wget --no-check-certificate https://osf.io/9xmh4/download -O subset_data.tar.gz
    tar -xzvf subset_data.tar.gz
    mkdir ill
    mv ./subset/subset_ill/* ill/
    mkdir ont/
    mv ./subset/subset_ont/* ont/
    mkdir rna/
    mv ./subset/subset_rna/* rna/
    """
}