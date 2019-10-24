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
                    path_db=\$(dirname !{db});
                    mkdir \$path_db/db/;
                    tar -xvf !{db} -C \$path_db/db/;
                    checkm data setRoot \$path_db/db;   
                    echo \$path_db/db/not_working > path_db.txt;
        fi
        """
    }