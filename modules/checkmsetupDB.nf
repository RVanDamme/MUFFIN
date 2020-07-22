process checkm_setup_db {
    label 'checkm'
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
                    echo -e "cat << EOF\\n\$path_db/db\\nEOF\\n" | checkm data setRoot   
                    echo \$path_db/db > path_db.txt;
        fi
        """
    }