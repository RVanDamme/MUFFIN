process checkm_setup_db {
    label 'checkm'
    errorStrategy = { task.exitStatus==14 ? 'retry' : 'terminate' }
    maxRetries = 5
    input:
    val(db)
    val(untar)
    output:
    path("path_db.txt")
    shell:
        """
        if [ !{untar} == true ] ;
        then
            checkm data setRoot !{db} ;
            echo '${db}' > path_db.txt;
        fi

        if [ !{untar} == false ] ;
        then
                    path_db=\$(dirname !{db});
                    mkdir -p \$path_db/db/;
                    tar -xvf ${db} -C \$path_db/db/;
                    checkm data setRoot \$path_db/db/ ;
                    echo \$path_db/db > path_db.txt;
        fi
        """
    }
