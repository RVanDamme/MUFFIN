process spades {
    label 'spades'

    //conda 'bioconda::spades=3.15.4'

    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/spades/", mode: 'copy', pattern: "assembly.fasta" 
    input:
    tuple val(name), path(illumina), path(ont)
    output:
    tuple val(name), path("assembly.fasta")
    
    script:
    """
    mem=\$(echo ${task.memory.toGiga()} | awk '{print int(\$1)}')
    cpus=\$(echo ${task.cpus})
    echo \$cpus \$mem
    spades.py -1 ${illumina[0]} -2 ${illumina[1]}  --meta --nanopore ${ont} -o spades_output -t \$cpus -m \$mem
    mv spades_output/contigs.fasta  assembly.fasta
    """

}

//mem=\$(echo ${task.memory} |sed 's/ GB//g' | sed 's/g//g')

process spades_short {
    label 'spades'

    conda 'bioconda::spades=3.15.4'
    
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/${name}/assemble/assembly/spades/", mode: 'copy', pattern: "assembly.fasta" 
    input:
    tuple val(name), path(illumina)
    output:
    tuple val(name), path("assembly.fasta")
    
    script:
    """
    mem=\$(echo ${task.memory} |sed 's/ GB//g' | sed 's/g//g')
    cpus=\$(echo ${task.cpus})
    echo \$cpus \$mem
    spades.py -1 ${illumina[0]} -2 ${illumina[1]} -o spades_output -t \$cpus -m \$mem
    mv spades_output/contigs.fasta  assembly.fasta
    """
}

// process spades {
//     label 'spades'
//     errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
//     maxRetries = 5
//     publishDir "${params.output}/${name}/assemble/assembly/spades/", mode: 'copy', pattern: "assembly.fasta" 
    
//     input:
//     tuple val(name), path(illumina), path(ont)
    
//     output:
//     tuple val(name), path("assembly.fasta")
    
//     script:
//     """
//     mem=\$(echo \$((${task.memory.toMega()} / 1024)) | awk '{print int(\$1)}')
//     cpus=\$(echo ${task.cpus})
//     echo "Using \$cpus CPUs and \$mem GB of memory for SPAdes."
    
//     # Define the base SPAdes command
//     spades_cmd="spades.py -1 ${illumina[0]} -2 ${illumina[1]} -o spades_output -t \$cpus -m \$mem"
    
//     # Add the nanopore reads to the command if they are provided
//     [[ -f "${ont}" ]] && spades_cmd="\$spades_cmd --nanopore ${ont}"
    
//     # Execute the SPAdes command
//     echo "Running: \$spades_cmd"
//     eval \$spades_cmd
    
//     mv spades_output/contigs.fasta assembly.fasta
//     """
// }

// if [ -n "${ont}" ]; then
//       spades_cmd="\${spades_cmd} --meta --nanopore ${ont}"
// fi