# MAFIN
Metagenomic Assembly pipeline using nextFlow for Illumina and Nanopore reads


## Introduction

MAFIN aims at being a reproducible pipeline for metagenome assembly
of crossed illumina and nanopore reads.

MAFIN uses the following software

| Task | Software | Version | Docker | integration|
| --- | --- | --- | --- | --- |
| QC illumina | [fastp](https://github.com/OpenGene/fastp) |  |  |  |
| QC ont | automated way to discard shortest reads |  |  |  |
|  | [filtlong](https://github.com/rrwick/Filtlong) |  |  |  |
| metagenomic composition of ont | [sourmash](https://sourmash.readthedocs.io/en/latest/) |  |  |  |
| Hybrid assembly | [Meta-spades](http://cab.spbu.ru/software/spades/) |  |  |  |
|  | [unicycler](https://github.com/rrwick/Unicycler) |  |  |  |
| Long read assembly | [MetaFlye](https://github.com/fenderglass/Flye) |  |  |  |
| polishing | [racon](https://github.com/isovic/racon) |  |  |  |
|  | [medaka](https://github.com/nanoporetech/medaka) |  |  |  |
|  | [pilon](https://github.com/broadinstitute/pilon/wiki) |  |  |  |
| mapping | [minimap2](https://github.com/lh3/minimap2) |  |  |  |
|  | [bwa](http://bio-bwa.sourceforge.net/) |  |  |  |
|  | [samtools](http://www.htslib.org/) |  |  |  |
| retrieve reads mapped to contig | [seqtk](https://github.com/lh3/seqtk) |  |  |  |
| Binning | [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) |  |  |  |
|  | [maxbin2](https://sourceforge.net/projects/maxbin2/) |  |  |  |
|  | [concoct](https://github.com/BinPro/CONCOCT) |  |  |  |
|  | [metawrap](https://github.com/bxlab/metaWRAP) |  |  |  |
| qc binning | [checkm](https://ecogenomics.github.io/CheckM/) |  |  |  |


## Installation

At the moment to install this pipeline and run this pipeline you need to use the conda installation:
```sh
#install the pipeline
git clone https://github.com/RVanDamme/MAFIN.git

#create an env and install metawrap
conda -y -p /path/to/install/metawrap-env python=2.7
source activate /path/to/install/metawrap-env
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels ursky
conda install -y -c ursky metawrap-mg
conda deactivate

#edit MAFIN/modules/metawrap_refine_bin.nf to use the env of metawrap
#you need to change the line 3 and 25 to the path of your env (/path/to/install/metawrap-env)
```

## Usage
the current 2 default usage are:
Spades
```
nextflow run MAFIN/main.nf --output results --assembler metaspades  --illumina fastq_ill/ --nanopore fastq_nano/ --core 2 --memory 16g  -profile conda
```
Flye
```
nextflow run MAFIN/main.nf --output results --assembler metaflye  --illumina fastq_ill/ --nanopore fastq_nano/ --core 2 --memory 16g  -profile conda
```
A more advanced  usage
```
nextflow run MAFIN/main.nf --output results --assembler metaspades  --illumina fastq_ill/ --nanopore fastq_nano/ --core 24 --memory 80g --out_qc --out_metawrap --out_unmapped -profile conda 
```
### Options

#### --cpus
* number of thread available
* default 2

#### --mem
* number of memory available
* default 16g

#### --assembler
* which method to use
* can be Hybrid with metaspades or long read + polishing with metaflye
* required

#### --illumina
* location of the dir containing the forward and reverse illumina reads in fasta or fastq 
* required

#### --nanopore
* location of the dir containing the nanopore reads in fasta or fastq
* required

#### --output
* output directory
* required

### Complete help and options
```
    *********Metagenomic Assembly pipeline using nextFlow for Illumina and Nanopore reads*********

    Mafin is composed of 2 part the retrieval of potential genome and the analysis of said genomes

        Usage example for retrieval:
    nextflow run mafin --retrieve --ont /path/to/ont_dir --illumina /path/to/illumina_dir --metaspades -profile conda
    or 
    nextflow run mafin --retrieve --ont /path/to/ont_dir --illumina /path/to/illumina_dir --metaflye -profile conda

        Input:
    --ont                       path to the directory containing the nanopore read file (fastq)
    -- illumina                 path to the directory containing the illumina read file (fastq)

        Output (default output is reassemblies from each bins):
    --output                    path to the output directory (default: $params.output)
    --assembly                  output the original assembly contigs file (default: false)
    --out_qc                    output the reads file after qc (default: false)
    --out_metabat               output the bins produce by metabat2 (default: false)
    --out_concoct               output the bins produce by concoct (default: false)
    --out_maxbin                output the bins produce by meaxbin2 (default: false)
    --out_metawrap              output the bins produce by metawrap refining (default: false)
    --out_bin_reads             output fastq files containing the reads mapped to each bin (default: false)
    --out_unmapped              output sorted bam files containing the unmmaped reads of illumina and nanopore (default:false)


    

        Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB for --metamaps [default: $params.memory]
    
        Options:
    --checkm_db                 path to an already INSTALLED checkm database (not the tar file)
    --checkm_tar_db             path to the tar checkm database (it will extract it in the dir)
    --sourmash                  path to an already installed sourmash database
    --skip_ill_qc               skip quality control of illumina files
    --skip_ont_qc               skip quality control of nanopore file
    --short_qc                  minimum size of the reads to be kept (default: $params.short_qc )
    --filtlong                  use filtlong to improve the quality furthermore (default: false)
    --model                     the model medaka will use (default: r941_min_high)
    --polish_iteration          number of iteration of pilon in the polish step (advanced)
    --extra_ill                 a list of additional ill sample file (with full path with a * instead of _R1,2.fastq) to use for the binning in Metabat2 and concoct
    --extra_ont                 a list of additional ont sample file (with full path) to use for the binning in Metabat2 and concoct
    --SRA_ill                   a list of additional ill sample from SRA accession number to use for the binning in Metabat2 and concoct
    --SRA_ont                   a list of additional ont sample from SRA accession number to use for the binning in Metabat2 and concoct
    --skip_metabat2             skip the binning using metabat2 (advanced)
    --skip_maxbin2              skip the binning using maxbin2 (advanced)
    --skip_concoct              skip the binning using concoct (advanced)

        Nextflow options:
    -profile                    change the profile of nextflow (currently available conda)
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)
```

## License

Code is [GPL-3.0](LICENSE)

## Contributing

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines
