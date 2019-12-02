The dammit modules were used to annotate the rna transcript using multiple Database
    they are deprecated because of multiple issues:
        - recent version(1.1) are not working
        - the installation of the software through nextflow -> conda is not working on test server
        - the database cannot be installed easily to a specify dir and when it's done the database is not recognise by dammit
    Instead eggNOG is used until better version of dammit