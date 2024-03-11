process semibin2 {
    maxForks 1
    label 'semibin2'

    conda 'bioconda::semibin=2.0.2'

    publishDir "\${params.output}/\${name}/assemble/binning/semibin2/", mode: 'copy', pattern: "bins_dir"
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    tuple val(name), path(assembly), path(bam_files)
    output:
    tuple val(name), path("bins_dir")
    
    script:
    """
    # Définition des options en fonction du choix du modèle
    model_option=""
    if ("${params.model}" == "pre-trained") {
        model_option="--environment ${params.environment}"
        echo "Running SemiBin2 with the pre-trained model: ${params.environment}"
    } else if ("${params.model}" == "self-supervised") {
        model_option="--self-supervised"
        echo "Running SemiBin2 with the self-supervised model. This process may take some time and requires disk space."
    } else if ("${params.model}" == "semi-supervised") {
        model_option="--semi-supervised"
        echo "Warning: Running SemiBin2 with the semi-supervised model is not recommended due to the long execution time and high memory requirement (~40GB). Consider using the self-supervised model instead."
    } else {
        echo "Defaulting to the pre-trained global model for SemiBin2."
        model_option="--environment ${params.environment}"
    }

    SemiBin2 single_easy_bin -i \${assembly} -b *.bam -o bins_dir/semibin2_bins -t \${task.cpus} \${model_option}
    """
}
