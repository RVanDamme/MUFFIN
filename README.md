# MAFIN
Metagenomic Assembly pipeline using nextFlow for Illumina and Nanopore reads


## Introduction

MAFIN aims at being a reproducible pipeline for metagenome assembly
of crossed illumina and nanopore reads.

MAFIN uses the following software

| Task | Software | Version | Docker | integration|
| --- | --- | --- | --- | --- |
| Filter of unwanted genome |  |  |  |  |
| assembly | [Meta-spades](http://cab.spbu.ru/software/spades/) |  |  |  |
| | [MetaFlye](https://github.com/fenderglass/Flye)/[Polishing] |  | | |
| QC assembly | [metaQUAST?](http://bioinf.spbau.ru/metaquast)| |   |  |
| mapping assembly | [???] |  |  | |
| Binning | [Metabat](https://bitbucket.org/berkeleylab/metabat/src/master/) |  |  | |

## Installation


## Usage


### Options

#### --cpus
* number of thread available
* default 4

#### --mem
* number of memory available
* default 8GB

#### --assembler
* which method to use
* can be Hybrid with metaspades or long read + polishing with metaflye
* required

#### --illumina_1 & --illumina_2
* location of the forward and reverse illumina reads in fasta or fastq 
* required

#### --nanopore
* location of the nanopore reads in fasta or fastq
* required

#### --output
* output directory
* required


## License

Code is [GPL-3.0](LICENSE)

## Contributing

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines
