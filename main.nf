#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
*********Start running MUFFIN*********
MUFFIN is a hybrid assembly and differential binning workflow for metagenomics, transcriptomics and pathway analysis.

If you use MUFFIN for your research please cite:

https://www.biorxiv.org/content/10.1101/2020.02.08.939843v1 

or

Van Damme R., Hölzer M., Viehweger A., Müller B., Bongcam-Rudloff E., Brandt C., 2020
"Metagenomics workflow for hybrid assembly, differential coverage binning, transcriptomics and pathway analysis (MUFFIN)",
doi: https://doi.org/10.1101/2020.02.08.939843 
**************************************
"""


XX = "21"
YY = "04"
ZZ = "0"

if ( nextflow.version.toString().tokenize('.')[0].toInteger() < XX.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}
else if ( nextflow.version.toString().tokenize('.')[1].toInteger() ==
XX.toInteger() && nextflow.version.toString().tokenize('.')[1].toInteger() <
YY.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}



include { helpMSG } from './modules/help.nf'
include { long_read_workflow } from './readflows/long.nf'
include { short_read_workflow } from './readflows/short.nf'
include { hybrid_workflow } from './readflows/hybrid.nf'

// Early exit and help message handling
if (params.help) { exit 0, helpMSG() }

// Error handling for execution profiles
if (!workflow.profile.matches(/^(?:local|sge|slurm|gcloud|ebi|lsf|git_action),(?:singularity|docker|conda)$/)) {
    exit 1, "Invalid execution profile. Please use -profile with a valid executer and engine combination, e.g., 'local,docker'."
}

// Always include readme and test modules
include {readme_output} from './modules/readme_output' params(output: params.output)
include {test} from './modules/test_data_dll'

//chackm2 database idea
params.db_path = '/home/user/databases'
params.db_file = 'uniref100.KO.1.dmnd'

workflow {

    // Validate assembler parameter
    if (!['metaflye', 'metaspades'].contains(params.assembler) &&  params.mode == "hybrid") {
        error "--assembler: ${params.assembler} is invalid. Should be 'metaflye' or 'metaspades' in hybrid mode."
    }
    
    switch (params.mode) {
        case "short":
            short_read_workflow()
            break
        case "long":
            long_read_workflow()
            break
        case "hybrid":
            hybrid_workflow()
            break
        default:
            error """Invalid analysis mode: ${params.mode}. Please use -mode with a valid mode: 
            "--mode short" for Illumina reads 
            "--mode long" for Nanopore reads 
            "--mode hybrid" for both Illumina and Nanopore"""
    }

}

