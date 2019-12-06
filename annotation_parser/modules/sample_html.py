#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
from urllib.request import urlopen
from urllib.error import HTTPError

def write_html_sample(dict_global_sample, output,
             globalpathwaylist, binnamelist, rna_pathway_list,
             dictrna):

    out = output+"/MAFIN_sample_result.html"
    outfile = open(out, "w")

    outfile.write("""
    <!doctype html>
    <html lang="en-US">
    <head>
        <meta charset="utf-8">
        <title>MAFIN Sample result</title>
        <meta name="author" content="Renaud Van Damme">
    </head>"""
    )

    num_bins = len(binnamelist)
    num_path = len(globalpathwaylist)
    num_path_rna = len(rna_pathway_list)
    outfile.write(f"""
    <body>
    <div id='summary'>
    <h1>Summary</h1>
    <h2>
    <ul>
        <li>Total number of bins: {num_bins}</li>
        <li>Total number of unique pathways in bins: {num_path}</li>
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
                    <th class="header">Bins [number of gene, number of expressed gene]</th>
                    <th class="header">Kegg Modules</th>
                </tr>
    """)

    for pathway in dict_global_sample:
        request = 'http://rest.kegg.jp/get/'+pathway
        print(request)
        req = urlopen(request)
        data = req.read().decode()
        data = data.split('\n')[1]
        pathway_name = re.split("NAME\s+", data)[1]
        list_gene= set()
        try:
            for gene in dictrna[pathway]:
                list_gene.add(gene)
                list_active_gene = "+".join(list(list_gene))
        except KeyError:
            list_active_gene=""

        outfile.write(f"""
        <tr>
        <td class="pathway_gene"><a href="https://www.kegg.jp/pathway/{pathway}+{list_active_gene}">{pathway_name}</a></td>
        <td class="pathway_gene">"""
        )
        for bins in dict_global_sample[pathway]:
            outfile.write(f"""{bins}[{dict_global_sample[pathway][bins][0]} 
                      ,{dict_global_sample[pathway][bins][2]}]; """)
        outfile.write("""</td>
        <td class="modules">
        """)
        mods=set()
        for bins in dict_global_sample[pathway]:
            for modules in dict_global_sample[pathway][bins][4]:
                mods.add(modules)
            for mod in mods:
                if mod != "":
                    outfile.write(f"""<a href='https://www.kegg.jp/module/{mod}+{list_active_gene}'>
                    {mod}</a>; 
                    """)

        outfile.write("""</td>
        </tr>
        """
        )
    outfile.close()
