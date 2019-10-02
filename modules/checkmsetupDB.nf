process checkm_download_db {
        label 'checkm'
        label 'ubuntu' 
    input:
    val(db)
    val(untar)
    shell:
        """
        if [ !{untar} == true ] ;
        then
            checkm data setRoot !{db} ;
        fi

        if [ !{untar} == false] ;
        then
            if test -f nextflow-autodownload-databases/checkm/checkm_data_2015_01_16.tar.gz ;
                then
                    tar -xvf nextflow-autodownload-databases/checkm/checkm_data_2015_01_16.tar.gz ;
                    checkm data setRoot nextflow-autodownload-databases/checkm/ ;
                
                else 
                    path_db=\$(dirname !{db});
                    tar -xvf !{db};
                    checkm data setRoot \$path_db;
            fi;
        fi
        """
    }