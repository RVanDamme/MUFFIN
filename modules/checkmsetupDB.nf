process checkm_setup_db {
    label 'checkm'
    input:
    val(db)
    val(untar)
    output:
    file("path_db.txt")
    shell:
        """
        if [ !{untar} == true ] ;
        then
            checkm data setRoot !{db} ;
            echo '!{db}' > path_db.txt;
        fi

        if [ !{untar} == false ] ;
        then
            if test -f nextflow-autodownload-databases/checkm/checkm_data_2015_01_16.tar.gz ;
                then
                    if [ -d "nextflow-autodownload-databases/checkm/db" ]  ;
                        then
                        if [ -L "nextflow-autodownload-databases/checkm/db" ]  ;
                            then
                                tar -xvf nextflow-autodownload-databases/checkm/checkm_data_2015_01_16.tar.gz -C nextflow-autodownload-databases/checkm/db_cheval ;
                                checkm data setRoot nextflow-autodownload-databases/checkm/db ;
                                echo "nextflow-autodownload-databases/checkm/db" > path_da.txt;
                            else
                                mkdir nextflow-autodownload-databases/checkm/db_poney/;
                                tar -xvf nextflow-autodownload-databases/checkm/checkm_data_2015_01_16.tar.gz -C nextflow-autodownload-databases/checkm/db_poney ;
                                checkm data setRoot nextflow-autodownload-databases/checkm/db ;
                                echo "nextflow-autodownload-databases/checkm/db" > path_da.txt;
                        fi;
                    fi;
                
                else 
                    path_db=\$(dirname !{db});
                    mkdir \$path_db/db/;
                    tar -xvf !{db} -C \$path_db/db/;
                    checkm data setRoot \$path_db/db;   
                    echo \$path_db/db > path_db.txt;
            fi;
        fi
        """
    }