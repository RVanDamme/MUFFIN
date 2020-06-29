process readme_output {
    label 'ubuntu'
    publishDir "${params.output}/", mode: 'copy', pattern: "README_output.txt"
    input:
    output:
        file("README.txt") 
    script:
        """
        wget https://osf.io/a6hru/download -O README_output.txt
        """
    }
}