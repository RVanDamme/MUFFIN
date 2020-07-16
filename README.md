# MUFFIN
MUFFIN is a hybrid assembly and differential binning workflow for metagenomics, transcriptomics and pathway analysis.

## INDEX

1. [Introduction](#introduction)
2. [Figure](#figure) :
    - [Workflow](#the-workflow)
    - [Parser output](#the-parser-output)
3. [Installation](#installation) :
    - [base installation](#base-installation)
    - [conda usage](#for-conda-usage)
    - [containers usage](#for-containers-usage-(docker/singularity))
    - [software installe locally](#for-usage-of-software-installed-locally)
4. [Usage](#usage) :
    - [Basic usage](#basic-usage)
    - [Advanced usage](#advanced-usage)
    - [Test the pipeline](#test-the-pipeline)
5. [Troubleshooting](#troubleshooting)
6. [Options](#options)
7. [Complete help and options](#complete-help-and-options)
8. [Bibliography](#bibliography)
9. [License](#license)

## Introduction

MUFFIN aims at being a reproducible pipeline for metagenome assembly
of crossed illumina and nanopore reads.

MUFFIN uses the following software

| Task | Software | Version | Docker | Image version|
| --- | --- | --- | --- | --- |
| QC illumina | [fastp](https://github.com/OpenGene/fastp) | 0.20.0 | [LINK](https://hub.docker.com/r/nanozoo/fastp) | 0.20.0--78a7c63 |
| QC ont | automated way to discard shortest reads |  |  |  |
|  | [filtlong](https://github.com/rrwick/Filtlong) | 0.2.0 | [LINK](https://hub.docker.com/r/nanozoo/filtlong) | v0.2.0--afa175e |
| metagenomic composition of ont | [sourmash](https://sourmash.readthedocs.io/en/latest/) | 2.0.1 | [LINK](https://hub.docker.com/r/nanozoo/sourmash) | 2.0.1--6970ddc |
| Hybrid assembly | [Meta-spades](http://cab.spbu.ru/software/spades/) | 3.13.1 | [LINK](https://hub.docker.com/r/nanozoo/spades) | 3.13.1--2c2a4c0 |
|  | [unicycler](https://github.com/rrwick/Unicycler) | 0.4.7 | [LINK](https://hub.docker.com/r/nanozoo/unicycler) | 0.4.7-0--c0404e6 |
| Long read assembly | [MetaFlye](https://github.com/fenderglass/Flye) | 2.7 | [LINK](https://hub.docker.com/r/nanozoo/flye) | 2.7--957a1a1 |
| polishing | [racon](https://github.com/lbcb-sci/racon) | 1.4.13 | [LINK](https://hub.docker.com/r/nanozoo/racon) | 1.4.13--bb8a908 |
|  | [medaka](https://github.com/nanoporetech/medaka) | 1.0.3 | [LINK](https://hub.docker.com/r/nanozoo/medaka) | 1.0.3--7c62d67 |
|  | [pilon](https://github.com/broadinstitute/pilon/wiki) | 1.23 | [LINK](https://hub.docker.com/r/nanozoo/pilon) | 1.23--b21026d |
| mapping | [minimap2](https://github.com/lh3/minimap2) | 2.17 | [LINK](https://hub.docker.com/r/nanozoo/minimap2) | 2.17--caba7af |
|  | [bwa](http://bio-bwa.sourceforge.net/) | 0.7.17 | [LINK](https://hub.docker.com/r/nanozoo/pilon) | 1.23--b21026d |
|  | [samtools](http://www.htslib.org/) | 1.9 | [LINK](https://hub.docker.com/r/nanozoo/minimap2) | 2.17--caba7af |
| retrieve reads mapped to contig | [seqtk](https://github.com/lh3/seqtk) | 1.3 | [LINK](https://hub.docker.com/r/nanozoo/seqtk) | 1.3--dc0d16b |
| Binning | [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) | 2.13 | [LINK](https://hub.docker.com/r/nanozoo/metabat2) | 2.13--0e2577e |
|  | [maxbin2](https://sourceforge.net/projects/maxbin2/) | 2.2.7 | [LINK](https://hub.docker.com/r/nanozoo/maxbin2) | 2.2.7--b643a6b |
|  | [concoct](https://github.com/BinPro/CONCOCT) | 1.1.0 | [LINK](https://hub.docker.com/r/nanozoo/concoct) | 1.1.0--03a3888 |
|  | [metawrap](https://github.com/bxlab/metaWRAP) | 1.2.2 | [LINK](https://hub.docker.com/r/nanozoo/metawrap) | 1.2.2--de94241 |
| qc binning | [checkm](https://ecogenomics.github.io/CheckM/) | 1.0.13 | [LINK](https://hub.docker.com/r/nanozoo/nanoplot) | 1.0.13--248242f |
|Taxonomic Classification  | [sourmash](https://sourmash.readthedocs.io/en/latest/) using the [gt-DataBase](https://gtdb.ecogenomic.org/) | 2.0.1 | [LINK](https://hub.docker.com/r/nanozoo/sourmash) | 2.0.1--6970ddc |
|  | [GTDB](https://gtdb.ecogenomic.org/) | version r89 |  |  |
| Annotations (bin and RNA) | [eggNOG](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2) | 2.0.1 | [LINK](https://hub.docker.com/r/nanozoo/eggnog-mapper) | 2.0.1--d5e0c8c |
|  | [eggNOG DB](http://eggnog5.embl.de/#/app/home) | v5.0 |  |  |
| *De novo* transcript and quantification | [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) | 2.9.1 | [LINK](https://hub.docker.com/r/nanozoo/trinity) | 2.9.1--82fe26c |
|  | [Salmon](https://github.com/COMBINE-lab/salmon) | 0.15.0 | [LINK](https://hub.docker.com/r/nanozoo/trinity) | 2.9.1--82fe26c |

## Figure

### The Workflow


![MUFFIN FLOWCHART FIGURE](.figure/Muffin_Workflow_simple.png)

### The parser output

![PARSER OUTPUT FIGURE](.figure/PANKEGG_simple.png)



## Installation

### base installation
You need to install nextflow Version 20.01+ ( https://www.nextflow.io/ )
```sh
# verify Java version (at least version 8+)
java -version 

# Setup nextflow (it will create a nextflow executable file in the current directory)
curl -s https://get.nextflow.io | bash

# If you want the pipeline installed locally use the following
git clone https://github.com/RVanDamme/MAFIN.git

# If you want to not install the pipeline use the following when running nextflow

nextflow run  RVanDamme/MUFFIN --parameters.....

```

### For conda usage
If you use conda you need to install Metawrap in an environment you create yourself, this is due to a known issue that will be fixed soon.

```sh

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

### For containers usage (docker/singularity)
If you use containers, you don't need extra installations

### For usage of software installed locally
You just need to have all the software used in the pipeline (see table above) installed and in your $PATH


## Usage

### Basic usage

```
path/to/nextflow run $MUFFIN_pipeline --output results_dir --assembler $assembler --illumina fastq_ill/ --ont fastq_ont/ --cpus 16 --memory 64g --modular full -profile $profile_executor,$profile_engine
```
 $MUFFIN_pipeline is either "path/to/MUFFIN/main.nf" or "RVanDamme/MUFFIN"

 $assembler is either:
  - "metaspades" for hybrid assembly
  - "metaflye" for long-read assembly with short-reads polishing

 $profile_executor can be:
  - "local" to run on your computer
  - "gcloud" to run on google life science cloud computing (you need to setup your project in nextflow.config)
  - "slurm" to run on HPC using slurm (e.g. UPPMAX)

 $profile_engine can be:
  - "local_engine" to execute the software installed locally
  - "conda" to execute using conda installation
  - "docker" and "singularity" to execute using the docker container

You can add "-resume" at the end of the command to restart it while keeping the process that succeeded
This is often used in case or error in the pipeline or if you modify slightly the command and want to avoid running everything again.
One exemple is to run the pipeline without RNA and rerun it adding RNA data, in this specific case the second time you add:
```
--rna path/to/rna -resume
```
Only the transcript processes and final parsing will be run

### Advanced usage

You can use RNA data to have transcriptomics analysis to do so add "--rna path/to/fastq_rna/"

You can run only partially the pipeline to do so change --modular to the right parameter:
 - "full" run the 3 steps of MUFFIN (assemble, classify, annotate)
 - "assemble" run the assembly and binning
 - "classify" run the classification of the bins (require a different input)
 - "annotate" run the annotation step of the bins (require a different input) and RNA if provided
 - "assemb-class" run assemble and classify step
 - "assemb-annot" run assemble and classify step
 - "class-annot" run classify and annotate step (require a different input)

To run classify and annotate independently from the assemble you need to provide a CSV file of the bins
The structure of the file should correspond to:
```
Samplename,path/to/bin1.fa
Samplename,path/to/bin2.fa
...
Samplename,path/to/binX.fa
Samplename,path/to/binY.fa
```
If you run "classify" with or without "annotation" use "--bin_classify"
If you run "annotate" without "classify" use "--bin_annotate"

### Test the pipeline
To test the pipeline we have a subset of 5 bins available at https://osf.io/9xmh4/
To run it you just need to add "test" in the -profile e.g.:
```
nextflow run RVanDamme/MUFFIN --output results_dir --assembler metaspades --cpus 8 --memory 32g --modular full -profile gcloud,docker,test
```
The subset contains also RNA data to test with transcriptomics analysis you just need to activate it using "--rna"
The results of the different test run are available at https://osf.io/m5czv/

## Troubleshooting
* If metawrap fail using conda check that you installed metawrap in a conda environment and put the path in "modules/metawrap_refine_bin.nf"

* If you run the pipeline with google life sciences and get error code 14
  It means the process was killed by google, you just need to run the pipeline again don't forget to add "-resume"

* If either Metawrap or checkm have the following error
  You need to increase the RAM in the command for local_engine and conda or in the "configs/containers.config"
```
IOError: [Errno 2] No such file or directory: 'binsA.checkm/storage/tree/concatenated.tre' 
``` 

* For other issue please open a ticket on ![github](https://github.com/RVanDamme/MUFFIN/issues)

## Options

### --cpus
* number of thread available
* default 2

### --mem
* number of memory available
* default 16g

### --assembler
* which method to use
* can be Hybrid with metaspades or long read + polishing with metaflye
* required

### --illumina
* location of the dir containing the forward and reverse illumina reads in fasta or fastq 
* required

### --nanopore
* location of the dir containing the nanopore reads in fasta or fastq
* required

### --output
* output directory
* required

### -profile
* engine and executor
* required
* engine: local_engine, conda, docker, singularity
* executor: local, gcloud, slurm


## Complete help and options
```
    *********hybrid assembly and differential binning workflow for metagenomics, transcriptomics and pathway analysis*********

    MUFFIN is still under development please wait until the first non edge version realease before using it.
    Please cite us using https://www.biorxiv.org/content/10.1101/2020.02.08.939843v1

    Mafin is composed of 3 part the assembly of potential metagenome assembled genomes (MAGs); the classification of the MAGs; and the annotation of the MAGs.

        Usage example:
    nextflow run RVanDamme/MUFFIN --output result --ont nanopore/ --illumina illumina/ --assembler metaspades --rna rna/ -profile local,docker
    or 
    nextflow run RVanDamme/MUFFIN --output result --ont nanopore/ --illumina illumina/ --assembler metaflye -profile local,docker

        Input:
    --ont                       path to the directory containing the nanopore read file (fastq) (default: $params.ont)
    --illumina                  path to the directory containing the illumina read file (fastq) (default: $params.illumina)
    --rna                       path to the directory containing the RNA-seq read file (fastq) (default: none)
    --bin_classify              path to the directory containing the bins files to classify (default: none)
    --bin_annotate              path to the directory containing the bins files to annotate (default: none)
    --assembler                 the assembler to use in the assembly step (default: $params.assembler)

        Optional input:
    --check_db                  path to the checkm database
    --check_tar_db              path to the checkm database tar compressed
    --sourmash_db               path to the LCA database for sourmash (default: GTDB LCA formated)
    --eggnog_db                 path to the eggNOG database

        Output:
    --output                    path to the output directory (default: $params.output)

        Outputed files:
        You can see the output structure at https://osf.io/a6hru/
    QC                          The reads file after qc
    Assembly                    The assembly contigs file 
    Bins                        The bins produced by CONCOCT, MetaBAT2, MaxBin2 and MetaWRAP (the refining of bins)
    Mapped bin reads            The fastq files containing the reads mapped to each metawrap bin
    Unmapped bin reads          The fastq files containing the unmmaped reads of illumina and nanopore
    Reassembly                  The reassembly files of the bins (.fa and .gfa)
    Checkm                      Various file outputed by CheckM (summary, taxonomy, plots and output dir)
    Sourmash                    The classification done by sourmash
    Classify summary            The summary of the classification and quality control of the bins (csv file)
    RNA output                  The de novo assembled transcript and the quantification by Salmon
    Annotation                  The annotations files from eggNOG (tsv format)
    Parsed output               HTML files that summarize the annotations and show graphically the pathways


    

        Basic Parameter:
    --cpus                      max cores for local use [default: $params.cpus]
    --memory                    80% of available RAM in GB for --metamaps [default: $params.memory]


        Workflow Options:
    --skip_ill_qc               skip quality control of illumina files
    --skip_ont_qc               skip quality control of nanopore file
    --short_qc                  minimum size of the reads to be kept (default: $params.short_qc )
    --filtlong                  use filtlong to improve the quality furthermore (default: false)
    --model                     the model medaka will use (default: $params.model)
    --polish_iteration          number of iteration of pilon in the polish step (default: $params.polish_iteration)
    --extra_ill                 a list of additional ill sample file (with full path with a * instead of _R1,2.fastq) to use for the binning in Metabat2 and concoct
    --extra_ont                 a list of additional ont sample file (with full path) to use for the binning in Metabat2 and concoct
    --skip_metabat2             skip the binning using metabat2 (advanced)
    --skip_maxbin2              skip the binning using maxbin2 (advanced)
    --skip_concoct              skip the binning using concoct (advanced)

        Nextflow options:
    -profile                    change the profile of nextflow both the engine and executor more details on github README
    -resume                     resume the workflow where it stopped
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)
```
## BIBLIOGRAPHY

BWA: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168] 

CheckM: Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. 2015. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25: 1043–1055.

Concoct: Johannes Alneberg, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie Schirmer, Joshua Quick, Umer Z Ijaz, Leo Lahti, Nicholas J Loman, Anders F Andersson & Christopher Quince. 2014. Binning metagenomic contigs by coverage and composition. Nature Methods, doi: 10.1038/nmeth.3103 

Fastp: Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

Filtlong: https://github.com/rrwick/Filtlong

Flye: Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, "Assembly of Long Error-Prone Reads Using Repeat Graphs", Nature Biotechnology, 2019 doi:10.1038/s41587-019-0072-8

HMMER: http://hmmer.org/ 

Maxbin2: Wu YW, Tang YH, Tringe SG, Simmons BA, and Singer SW, "MaxBin: an automated binning method to recover individual genomes from metagenomes using an expectation-maximization algorithm", Microbiome, 2:26, 2014.

Medaka: https://github.com/nanoporetech/medaka

Metabat2: Kang DD, Froula J, Egan R, Wang Z. MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities. PeerJ 2015;3:e1165. doi:10.7717/peerj.1165

Metawrap: Uritskiy, G.V., DiRuggiero, J. and Taylor, J. (2018). MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome, 6(1). https://doi.org/10.1186/s40168-018-0541-1

Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

Pilon: Bruce J. Walker, Thomas Abeel, Terrance Shea, Margaret Priest, Amr Abouelliel, Sharadha Sakthikumar, Christina A. Cuomo, Qiandong Zeng, Jennifer Wortman, Sarah K. Young, Ashlee M. Earl (2014) Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLoS ONE 9(11): e112963. doi:10.1371/journal.pone.0112963

pplacer: Matsen FA, Kodner RB, Armbrust EV. 2010. pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: doi:10.1186/1471-2105-11-538.

prodigal: Hyatt D, Locascio PF, Hauser LJ, Uberbacher EC. 2012. Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28: 2223–2230.

Racon: Vaser R, Sovic I, Nagarajan N, Sikic M. 2017. Fast and accurate de novogenome assembly from long uncorrected reads. Genome Res 27:737–746.https://doi.org/10.1101/gr.214270.116

Samtools: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9 [19505943]

Seqtk: https://github.com/lh3/seqtk

Sourmash: Brown et al, (2016), sourmash: a library for MinHash sketching of DNA, Journal of Open Source Software, 1(5), 27, doi:10.21105/joss.00027

Spades:  Lapidus A., Antipov D., Bankevich A., Gurevich A., Korobeynikov A., Nurk S., Prjibelski A., Safonova Y., Vasilinetc I., Pevzner P. A. New Frontiers of Genome Assembly with SPAdes 3.0.	(poster), 2014 

Unicycler: Wick RR, Judd LM, Gorrie CL, Holt KE (2017) Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595. https://doi.org/10.1371/journal.pcbi.1005595

## License

Code is [GPL-3.0](LICENSE)

## Contributing

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines
