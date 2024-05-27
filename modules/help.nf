def helpMSG() {
    log.info """
    ********* Hybrid Assembly and Differential Binning Workflow for Metagenomics, Transcriptomics, and Pathway Analysis *********

    Please cite us using: https://www.biorxiv.org/content/10.1101/2020.02.08.939843v1

    MUFFIN comprises three main parts: the assembly of potential metagenome-assembled genomes (MAGs), the classification of the MAGs, and the annotation of the MAGs.

    Usage Examples:
        nextflow run main.nf --ont nanopore/ --illumina illumina/ --assembler metaspades --rna rna/ -profile local,docker
        nextflow run main.nf --ont nanopore/ --illumina illumina/ --assembler metaflye -profile local,conda

    Input:
    --ont                       Path to the directory containing the nanopore read files (fastq). (Default: $params.ont)
    --illumina                  Path to the directory containing the illumina read files (fastq). (Default: $params.illumina)
    --rna                       Path to the directory containing the RNA-seq read files (fastq). (Default: none)
    --bin_classify              Path to the directory containing the bin files to classify. (Default: none)
    --bin_annotate              Path to the directory containing the bin files to annotate. (Default: none)
    --assembler                 Assembler to use in the assembly step. (Default: $params.assembler)

    Optional Input:
    --checkm2db                 Path to the CheckM database.
    --checkm2db_force_update    Force CheckM2 to download the diamond database.
    --sourmash_db               Path to the LCA database for Sourmash. (Default: gtdb-rs214-k31.lca.json)
    --sourmash_db_full          Path to the full database for Sourmash. (Default: gtdb-rs214-k31)
    --sourmash_db_lineage       Path to the lineage database for Sourmash. (Default: gtdb-rs214.lineages.csv)
    --eggnog_db                 Path to the eggNOG database.

    Output:
    --output                    Path to the output directory. (Default: $params.output)

    Output Files:
    You can see the output structure at: https://osf.io/a6hru/
    QC                          Reads files after quality control.
    Assembly                    Assembly contigs file.
    Bins                        Bins produced by MetaBAT2, Semibin2, and ComeBin.
    CheckM2                     Various files output by CheckM2 (summary, taxonomy, plots, and output directory).
    Unmapped bin reads          Fastq files containing the unmapped reads of illumina and nanopore.
    Sourmash                    Classification results from Sourmash.
    Classify summary            Summary of the classification and quality control of the bins (CSV file).
    RNA output                  De novo assembled transcript and quantification by Salmon.
    Annotation                  Annotation files from eggNOG (TSV format).

    Basic Parameters:
    --cpus                      Maximum cores for local use. (Default: $params.cpus)
    --memory                    Maximum available RAM in GB for --metamaps. (Default: $params.memory)

    Workflow Options:
    --skip_ill_qc               Skip quality control of illumina files.
    --skip_ont_qc               Skip quality control of nanopore files.
    --short_qc                  Minimum size of reads to be kept. (Default: $params.short_qc)
    --ont_min_qc                Minimum quality of reads to be kept. (Default: $params.ont_min_qc)
    --checkm2_low               Use low memory usage for CheckM2.
    --skip_bin_sorting          Skip bin sorting based on GTDB quality value "completeness - 5*contamination > 50".
    --skip_bad_reads_recovery   Skip retrieval and classification of unused reads after bin sorting.
    --skip_pilon                Skip polishing of bin sequences for long-read assemblies. (Default: $params.skip_pilon)
    --quast                     Enable MetaQUAST analysis. (Default: $params.quast)
    --model                     Model to be used by Medaka. (Default: $params.model)
    --polish_iteration          Number of iterations of Pilon in the polishing step. (Default: $params.polish_iteration)
    --bintool                   Select binning tool: metabat2, semibin2, or comebin. (Default: $params.bintool)

    Semibin2 Options (Advanced):
    --bining_model              Binning model for Semibin2. (Default: $params.bining_model)
                                Options: pre-trained, semi-supervised, self-supervised
    --environment               Select environment for Semibin2. (Default: $params.environment)
                                See: https://github.com/BigDataBiology/SemiBin

    Nextflow Options:
    -profile                    Change the profile of Nextflow for engine and executor settings. More details on GitHub README.
    -resume                     Resume the workflow from where it stopped.
    -with-report rep.html       Generate a report of CPU and RAM usage. (May cause errors)
    -with-dag chart.html        Generate a flowchart of the process tree.
    -with-timeline time.html    Generate a timeline of the workflow. (May cause errors)
    """
}
