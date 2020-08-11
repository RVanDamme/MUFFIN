process readme_output {
    label 'ubuntu'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    publishDir "${params.output}/", mode: 'copy', pattern: "README_output.txt"
    input:
    output:
        path("README_output.txt") 
    script:
        """
        wget https://osf.io/a6hru/download -O README_output.txt
        """

}