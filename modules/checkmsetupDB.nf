process checkm_download_db {
        label 'checkm'
        label 'ubuntu' 
    input:
    val(db)
    val(untar)
    shell:
        """
        if [ !{untar} == true ]
        then
            checkm data setRoot !{db}
        fi

        if [ !{untar} == false]
        then
            path_db=\$(dirname !{db})
            cd \$path_db
            tar -xvf !{db}
            checkm data setRoot \$path_db
        fi
        """
    }