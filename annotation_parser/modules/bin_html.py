import re
from urllib.request import urlopen
from urllib.error import HTTPError


def write_html_bins(dict_global_bin, output,
                      globalpathwaylist, rna_pathway_list,
                      dictrna):
    for bin_html in dict_global_bin.keys():
        out = output+"/MAFIN_"+bin_html+"_result.html"
        outfile = open(out, "w")

        outfile.write(f"""
        <!doctype html>
        <html lang="en-US">
        <head>
            <meta charset="utf-8">
            <title>MAFIN {bin_html} result</title>
            <meta name="author" content="Renaud Van Damme">
        </head>"""
                    )

        num_path_bin = len(dict_global_bin[bin_html])
        num_path = len(globalpathwaylist)
        num_path_rna = len(rna_pathway_list)
        outfile.write(f"""
        <body>
        <div id='summary'>
        <h1>Summary</h1>
        <h2>
        <ul>
            <li>Total number of unique pathway in this bin: {num_path_bin}</li>
            <li>Total number of unique pathways in all bins: {num_path}</li>
            <li>Total number of unique pathways in RNA: {num_path_rna}</li>
            <li>This file contains only the eggNOG annotation that have a kegg pathway id, for further research please look at the annotations.tsv files</li>
            <li>This result file was produced by <a href="https://github.com/RVanDamme/MAFIN">MAFIN</a> </li>
        </ul>
        </h2>
        </div>
        """
                    )
        outfile.write("""
        <h1>Pathways</h1>
            <style type="text/css">
                .tg {
                    border-collapse: collapse;
                    border-spacing: 0;
                }
            
                .tg td {
                    font-family: Arial, sans-serif;
                    font-size: 14px;
                    padding: 10px 5px;
                    border-style: solid;
                    border-width: 1px;
                    overflow: hidden;
                    word-break: normal;
                    border-color: black;
                }
            
                .tg th {
                    font-family: Arial, sans-serif;
                    font-size: 14px;
                    font-weight: normal;
                    padding: 10px 5px;
                    border-style: solid;
                    border-width: 1px;
                    overflow: hidden;
                    word-break: normal;
                    border-color: black;
                }
            
                .tg .header {
                    font-weight: bold;
                    text-decoration: underline;
                    font-size: x-large;
                    font-family: Georgia, serif !important;
                    ;
                    background-color: #82aaea;
                    border-color: inherit;
                    text-align: center;
                    vertical-align: top
                }
            
                .tg .pathway_gene {
                    font-size: medium;
                    font-family: Tahoma, Geneva, sans-serif !important;
                    ;
                    border-color: inherit;
                    text-align: center;
                    vertical-align: middle
                }
                .tg .modules {
                    font-size: small;
                    font-family: Tahoma, Geneva, sans-serif !important;
                    ;
                    border-color: inherit;
                    text-align: left;
                    vertical-align: middle
                }
            
                @media screen and (max-width: 767px) {
                    .tg {
                        width: auto !important;
                    }
            
                    .tg col {
                        width: auto !important;
                    }
            
                    .tg-wrap {
                        overflow-x: auto;
                        -webkit-overflow-scrolling: touch;
                    }
                }
            </style>
            <div class="tg-wrap">
                <table class="tg">
                    <tr>
                        <th class="header">Pathways</th>
                        <th class="header">Expressed genes</th>
                        <th class="header">Kegg Modules</th>
                        <th class="header">Non expressed genes</th>
                    </tr>
        """)


        for pathway in dict_global_bin[bin_html]:
                request = 'http://rest.kegg.jp/get/'+pathway
                print(request)
                req = urlopen(request)
                data = req.read().decode()
                data = data.split('\n')[1]
                pathway_name = re.split("NAME\s+", data)[1]
                set_active_gene = set()
                if dict_global_bin[bin_html][pathway][3] != "":
                    for gene in dict_global_bin[bin_html][pathway][3]:
                        set_active_gene.add(gene)
                set_gene = set()
                for gene in dict_global_bin[bin_html][pathway][1]:
                    set_gene.add(gene)
                if dict_global_bin[bin_html][pathway][3] != "":
                    list_html_gene = "+".join(list(set_active_gene))
                    list_active_gene = list(set_active_gene)
                    list_inactive_gene = [activ for activ in list(set_gene) if activ not in list_active_gene]
                else:
                    list_html_gene = ""
                    list_active_gene = ""
                    list_inactive_gene = list(set_gene)
                outfile.write(f"""
                <tr>
                <td class="pathway_gene"><a href="https://www.kegg.jp/pathway/{pathway}+{list_html_gene}">{pathway_name}</a></td>
                <td class="pathway_gene">"""
                )
                for gene in list_active_gene:
                    request_gene = 'http://rest.kegg.jp/get/'+gene
                    print(request_gene)
                    req_gene = urlopen(request_gene)
                    data_gene = req_gene.read().decode()
                    data_gene = data_gene.split('\n')[1]
                    gene_name = re.split("NAME\s+", data_gene)[1]
                    outfile.write(
                        f"""<a href="https://www.kegg.jp/dbget-bin/www_bget?{gene}">[{gene_name}]</a>; """)
                outfile.write("""</td>
                <td class="modules">
                """)
                mods = set()
                for modules in dict_global_bin[bin_html][pathway][4]:
                    mods.add(modules)
                for mod in mods:
                    if mod != "":
                        outfile.write(f"""<a href='https://www.kegg.jp/module/{mod}+{list_html_gene}'>
                        {mod}</a>;
                        """)
                outfile.write("""</td>
                    <td class="pathway_gene">
                    """)
                for gene in list_inactive_gene:
                    request_gene = 'http://rest.kegg.jp/get/'+gene
                    print(request_gene)
                    req_gene = urlopen(request_gene)
                    data_gene = req_gene.read().decode()
                    data_gene = data_gene.split('\n')[1]
                    gene_name = re.split("NAME\s+", data_gene)[1]
                    outfile.write(
                        f"""<a href="https://www.kegg.jp/dbget-bin/www_bget?{gene}">[{gene_name}]</a>; """)
                outfile.write("""</td>
                </tr>
                """
                )
        outfile.close()
