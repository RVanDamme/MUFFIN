#!/usr/bin/env python

import os
import re
import csv
import sys
import shutil
import logging
import argparse
from urllib.error import HTTPError 
# used to request the name of the pathways and genes
from urllib.request import urlopen
from collections import defaultdict


def rna_parse(rna, dict_transcript, dictrna, rna_pathway_list):
        pathway_list = []
        with open(rna) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        line_count += 1
                    else:
                        for gene in row[8].split(','):
                            if gene != '':
                                if not row[0] in dict_transcript:
                                        dict_transcript[row[0]] = [gene]
                                else:
                                        dict_transcript[row[0]].append(gene)
                        for pathway in row[9].split(','):
                            if pathway.startswith("map"):
                                pass
                            elif pathway.startswith("ko"):
                                pathway_list.append(pathway)
                                for gene in row[8].split(','):
                                    if gene != '':
                                        if not pathway in dictrna:
                                                dictrna[pathway] = [gene]
                                        else:
                                                dictrna[pathway].append(gene)
                            elif pathway == "":
                                pass
                            else:
                                linenumber = line_count+1
                                print(
                                    f"error unknown pathway {pathway} in line number {linenumber} in {rna}")
                        line_count += 1
        rna_pathway_list = list(set(pathway_list))
        return dict_transcript, dictrna, rna_pathway_list


def quant_parse(dictquant, level):
    with open(level) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                        line_count += 1
            else:
                if not row[0] in dictquant:
                    dictquant[row[0]] = [row[4]]
                else:
                    dictquant[row[0]].append(row[4])
                line_count += 1
        return dictquant


def genelevel(output, dict_transcript, dictquant):
    out = output+"/genelevel.csv"
    with open(out, mode='w', newline='') as gene_file:
        gene_write = csv.writer(gene_file, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        gene_write.writerow(["transcript id", "pathways", "Number of reads"])
        for transcript in dict_transcript.keys():
            try:
                gene_write.writerow(
                    [transcript, dict_transcript[transcript], dictquant[transcript]]
                )
            except KeyError:
                pass
    return


def bin_parse(bins,
              globalpathwaylist, binnamelist,
              dict_global_sample, dict_global_bin, dictrna):
    pathwaylist = []
    for binfile in bins:
        bin_name = os.path.splitext(os.path.splitext(
            os.path.basename(binfile))[0])[0]
        binnamelist.append(bin_name)

        listko = []
        listactiv = []
        dictko = {}
        dictactiv = {}
        dictmodules = {}

        with open(binfile) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                line_count = 0
                for row in csv_reader:
                    if line_count == 0:
                        line_count += 1
                    else:
                        for pathway in row[9].split(','):
                            if pathway.startswith("map"):
                                pass
                            elif pathway.startswith("ko"):
                                pathwaylist.append(pathway)
                                for gene in row[8].split(','):
                                    if gene != '':
                                        if not pathway in dictko:
                                            dictko[pathway] = [gene]
                                        else:
                                            dictko[pathway].append(gene)
                                        if pathway in dictrna:
                                            if gene in dictrna[pathway]:
                                                if not pathway in dictactiv:
                                                    dictactiv[pathway] = [gene]
                                                else:
                                                    dictactiv[pathway].append(
                                                        gene)
                                for module_number in row[10].split(','):
                                    if module_number != '':
                                        if not pathway in dictmodules:
                                            dictmodules[pathway] = [
                                                module_number]
                                        else:
                                            dictmodules[pathway].append(
                                                module_number)

                            elif pathway == "":
                                pass
                            else:
                                linenumber = line_count+1
                                print(
                                    f"error unknown pathway {pathway} in line number {linenumber} in {binfile}")
                        line_count += 1
        csv_file.close()
        for pathw in dictko.keys():
            number_ko = len(dictko[pathw])
            try:
                number_activ = len(dictactiv[pathw])
            except KeyError:
                number_activ = 0
            try:
                if number_activ != 0:
                    dict_global_sample[pathw][bin_name] = [
                        number_ko, dictko[pathw], number_activ, dictactiv[pathw], dictmodules[pathw]]
                    dict_global_bin[bin_name][pathw] = [
                        number_ko, dictko[pathw], number_activ, dictactiv[pathw], dictmodules[pathw]]
                else:
                    dict_global_sample[pathw][bin_name] = [
                        number_ko, dictko[pathw], number_activ, "", dictmodules[pathw]]
                    dict_global_bin[bin_name][pathw] = [
                        number_ko, dictko[pathw], number_activ, "", dictmodules[pathw]]
            except KeyError:
                if number_activ != 0:
                    dict_global_sample[pathw][bin_name] = [
                        number_ko, dictko[pathw], number_activ, dictactiv[pathw], ""]
                    dict_global_bin[bin_name][pathw] = [
                        number_ko, dictko[pathw], number_activ, dictactiv[pathw], ""]
                else:
                    dict_global_sample[pathw][bin_name] = [
                        number_ko, dictko[pathw], number_activ, "", ""]
                    dict_global_bin[bin_name][pathw] = [
                        number_ko, dictko[pathw], number_activ, "", ""]
    globalpathwaylist = list(set(pathwaylist))
    return dict_global_sample, dict_global_bin, globalpathwaylist, binnamelist


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
        list_gene = set()
        try:
            for gene in dictrna[pathway]:
                list_gene.add(gene)
                list_active_gene = "+".join(list(list_gene))
        except KeyError:
            list_active_gene = ""

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
        mods = set()
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
                    list_inactive_gene = [activ for activ in list(
                        set_gene) if activ not in list_active_gene]
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


def parse(args):
    """
    Main function of the parser

    Arguments: args (object): all the dictionnary argument from argparse
    """
    logger = logging.getLogger(__name__)
    output = args.output
    bins = args.bins[0].split()
    level = args.level
    rna = args.rna
    # try:
    os.makedirs(args.output)

    dict_transcript = {}
    dictrna = {}
    rna_pathway_list = []
    dict_transcript, dictrna, rna_pathway_list = rna_parse(
        rna, dict_transcript, dictrna, rna_pathway_list)

    dictquant = {}
    dictquant = quant_parse(dictquant, level)

    genelevel(output, dict_transcript, dictquant)

    globalpathwaylist = []
    binnamelist = []
    dict_global_sample = defaultdict(dict)
    dict_global_bin = defaultdict(dict)
    dict_global_sample, dict_global_bin, globalpathwaylist, binnamelist = bin_parse(bins,
                                                                                    globalpathwaylist, binnamelist,
                                                                                    dict_global_sample, dict_global_bin, dictrna)
    write_html_sample(dict_global_sample, output,
                      globalpathwaylist, binnamelist, rna_pathway_list,
                      dictrna)
    write_html_bins(dict_global_bin, output,
                    globalpathwaylist, rna_pathway_list,
                    dictrna)

    # except OSError as e:
    #     logger.error(f"{args.output} output directory already exists. Aborting.")
    #     sys.exit(1)
    # else:
    #     pass


def main():
    # argument handling
    parser = argparse.ArgumentParser(
        description='Parse annotations files to output nice HTML results')
    parser_input = parser.add_argument_group("required input")
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print version and exit"
    )
    if not '-v' in sys.argv:
        if not '--version' in sys.argv:
            parser_input.add_argument(
                "-b",
                "--bins",
                required=True,
                dest="bins",
                help="List of bins annotated by eggnog (required)",
                nargs='+',
                metavar='bin.?.annotation.tsv'
            )
            parser_input.add_argument(
                "-r",
                "--rna",
                required=True,
                dest="rna",
                help="Annotation file of the rna transcript (required)",
                metavar='sample.annotations.tsv'
            )
            parser_input.add_argument(
                "-l",
                "--level",
                required=True,
                dest="level",
                help="Quantification of the transcripts expression by trinity/salmon (required)",
                metavar='sample_transcript_quant.sf'
            )

    parser_output = parser.add_argument_group("output directory")
    parser_output.add_argument(
        "-o",
        "--output",
        nargs='?',
        const=1,
        default="./parser_results",
        dest="output", metavar="Out_dir_name",
        help="Name for the output directory"
    )

    args = parser.parse_args()

    try:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger("PANKEGG")
        if args.version:
            logger.info(f" version 1.0.0")
            sys.exit(0)
        parse(args)
    except AttributeError as e:
        logger.debug(e)
        parser.print_help()
        raise


if __name__ == "__main__":
    main()
