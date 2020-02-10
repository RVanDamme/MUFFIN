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


def bin_parse(bins,
              globalpathwaylist, binnamelist, dictgeneral,
              dict_global_sample, dict_global_bin):
    pathwaylist = []
    for binfile in bins:
        bin_name = os.path.splitext(os.path.splitext(
            os.path.basename(binfile))[0])[0]
        binnamelist.append(bin_name)

        listko = []
        dictko = {}
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
                                        if pathway not in dictko:
                                            dictko[pathway] = [gene]
                                        else:
                                            dictko[pathway].append(gene)
                                        if pathway not in dictgeneral:
                                            dictgeneral[pathway]= [gene]                                        
                                        elif pathway in dictgeneral:
                                            if gene not in dictgeneral[pathway]:
                                               dictgeneral[pathway].append(gene)
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
            number_general = len(dictgeneral[pathw])   
            dict_global_sample[pathw][bin_name] = [
                number_ko, dictko[pathw]]
            dict_global_bin[bin_name][pathw] = [
                number_general, number_ko, dictko[pathw]]
            dictgeneral[pathw]=list(set(dictgeneral[pathw]))   
    globalpathwaylist = list(set(pathwaylist))
    return dictgeneral, dict_global_sample, dict_global_bin, globalpathwaylist, binnamelist


def write_html_sample(dictgeneral, dict_global_sample, output,
                      globalpathwaylist, binnamelist):

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
    outfile.write(f"""
    <body>
    <div id='summary'>
    <h1>Summary</h1>
    <h2>
    <ul>
        <li>Total number of bins: {num_bins}</li>
        <li>Total number of unique pathways in bins: {num_path}</li>
        <li>This file contains only the eggNOG annotation that have a kegg pathway id, for further research please look at the annotations.tsv files</li>
        <li>This result file was produced by <a href="https://github.com/RVanDamme/MAFIN">MAFIN</a> </li>
    </ul>
    </h2>
    </div>
    """
    )
    outfile.write("""
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
              font-size: x-large;
              font-family: Georgia, serif !important;
              ;
              background-color: #7097AB;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .header2 {
              font-size: large;
              font-family: Georgia, serif !important;
              ;
              background-color: #7097AB;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .helpfont {
              font-size: medium;
              font-family: Georgia, serif !important;
              ;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .pathway_gene {
              font-size: medium;
              font-family: Tahoma, Geneva, sans-serif !important;
              ;
              background-color: #E8E8E8;
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
          .collapsible {
              background-color: #7097AB;
              color: black;
              cursor: pointer;
              padding: 18px;
              width: 100%;
              border: none;
              text-align: left;
              outline: none;
              font-size: large;
            }

        .active, .collapsible:hover {
              background-color: #555;
            }

        .content {
          padding: 0 18px;
          display: none;
          overflow: hidden;
          background-color: #f1f1f1;
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
    """)
    outfile.write(f"""
    <body>
        <div id='summary'>
        <h1>Help</h1>
        <button type="button" class="collapsible">Open Help</button>
        <div class="content"><p>
        <ol>
            <li>INDEX of the table
            <ul>
                <li>The columns "Pathways" is composed of name of each pathway present in the bins and link to a figure with the orthologs present in the bins in <font color="green">green</font></li>
                <li>The column "Bins Composition" is the list of the bins with orthologs present in the pathway plus for each bin the number of orthologs from the bin present in the pathway</li>
                <li>The column "Number of orthologs" is the total number of orthologs of the pathway present across all bins</li>
                <li>The column "Orthologs" is the list of all the orthologs of the pathway present across all bins</li>
            </ul></li>
            <li> Figure detail
            <ul>
                <li>The Figures in the links: <ul> 
                    <li>The orthologs present in the bins are in green</li>
            <li>Troubleshooting
            <ul>
                <li>When the link of the pathway is not loading or not showing anything, it means that there is too much orthologs to show on the figure.
            Try to strip everything after "https://www.kegg.jp/kegg-bin/show_pathway?PATWAY_ENTRY_NUMBER/" to still see the pathway</li>
            </ul></li>
        </ol>
        </p></div>
        </div>
    """
    )

    outfile.write("""
    <script>
        var coll = document.getElementsByClassName("collapsible");
        var i;

        for (i = 0; i < coll.length; i++) {
          coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            var content = this.nextElementSibling;
            if (content.style.display === "block") {
              content.style.display = "none";
            } else {
              content.style.display = "block";
            }
          });
        }
        </script>
    """)
    
    outfile.write("""
    <h1>Pathways</h1>
        
        <div class="tg-wrap">
            <table class="tg">
                <tr>
                    <th class="header">Pathways Summary</th>
                    <th class="header">Bins Compositions</th>
                    <th class="header">Number of orthologs</th>
                    <th class="header">Orthologs</th>
                </tr>
                <tr>
                    <th class="header2"><font color="#e7bcd4">▉▉</font>represent the orthologs present in the Bins</th>
                    <th class="header2">structure: Bins[<font color="#e7bcd4">▉▉</font>number of orthologs in the bin]</th>
                    <th class="header2">represent the total number of orthologs present in the pathway across the bins</th>
                    <th class="header2">list of the orthologs present in the bins</th>
                </tr>
    """)

    for pathway in dict_global_sample:
        request = 'http://rest.kegg.jp/get/'+pathway
        print(request)
        try:
            req = urlopen(request)
            data = req.read().decode()
            data = data.split('\n')[1]
            pathway_name = re.split("NAME\s+", data)[1]
            set_total_gene = set()
            try:
                for gene in dictgeneral[pathway]:
                    set_total_gene.add(gene)
                set_html= set()
                for gene in set_total_gene:
                    set_html.add(gene)
                list_html = "".join(set_html)
            except KeyError:
                list_html = ""
            outfile.write(f"""
            <tr>
            <td class="pathway_gene"><a href="https://www.kegg.jp/kegg-bin/show_pathway?{pathway}/{list_html}/default%3d%23e7bcd4">{pathway_name}</a></td>
            <td class="pathway_gene">"""
                        )
            
            for bins in dict_global_sample[pathway]:
                outfile.write(
                    f"""{bins}[<span style="background-color:#e7bcd4">{len(dict_global_sample[pathway][bins][1])}</span>]; """)
            n_gene=len(dictgeneral[pathway])
            outfile.write(f"""</td>
            <td class="pathway_gene">{n_gene}
            </td>
            <td class="listgenes">
            """)   
            for gene in dictgeneral[pathway]:
                request_gene = 'http://rest.kegg.jp/get/'+gene
                print(request_gene)
                try:
                    req_gene = urlopen(request_gene)
                    data_gene = req_gene.read().decode()
                    data_gene = data_gene.split('\n')[1]
                    gene_name = re.split("NAME\s+", data_gene)[1]
                    outfile.write(
                        f"""<a href="https://www.kegg.jp/dbget-bin/www_bget?{gene}">[{gene_name}]</a>; """)
                except:
                    outfile.write(f"""[{gene} unknow by KEGG DB]""")
            outfile.write("""</td>
            </tr>
            """
                        )
        except:
            outfile.write(f"""<tr>
			<td class="pathway_gene"> {pathway} unknow by the KEGG DATABASE </td></tr> """)
    outfile.close()


def write_html_bins(dictgeneral, dict_global_bin, output,
                    globalpathwaylist):
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
        outfile.write(f"""
        <body>
        <div id='summary'>
        <h1>Summary</h1>
        <h2>
        <ul>
            <li>Total number of unique pathway in this bin: {num_path_bin}</li>
            <li>Total number of unique pathways in all bins: {num_path}</li>
            <li>This file contains only the eggNOG annotation that have a kegg pathway id, for further research please look at the annotations.tsv files</li>
            <li>This result file was produced by <a href="https://github.com/RVanDamme/MAFIN">MAFIN</a> </li>
        </ul>
        </h2>
        </div>
        """
        )

        outfile.write("""
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
              font-size: x-large;
              font-family: Georgia, serif !important;
              ;
              background-color: #7097AB;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .header2 {
              font-size: large;
              font-family: Georgia, serif !important;
              ;
              background-color: #7097AB;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .helpfont {
              font-size: medium;
              font-family: Georgia, serif !important;
              ;
              border-color: inherit;
              text-align: center;
              vertical-align: top
          }
          .tg .pathway_gene {
              font-size: medium;
              font-family: Tahoma, Geneva, sans-serif !important;
              ;
              background-color: #E8E8E8;
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
          .collapsible {
              background-color: #7097AB;
              color: black;
              cursor: pointer;
              padding: 18px;
              width: 100%;
              border: none;
              text-align: left;
              outline: none;
              font-size: large;
            }

        .active, .collapsible:hover {
              background-color: #555;
            }

        .content {
          padding: 0 18px;
          display: none;
          overflow: hidden;
          background-color: #f1f1f1;
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
        """)
        outfile.write(f"""
        <body>
        <div id='summary'>
        <h1>Help</h1>
        <button type="button" class="collapsible">Open Help</button>
        <div class="content"><p>
        <ol>
            <li>INDEX of the table
            <ul>
                <li>The column "Pathways" is composed of name of each pathway present in the bins and link to a figure with the orthologs present in this bin in <font color="#e7bcd4">▉▉</font></li>
                <li>The column "Orthologs" is the list of the orthologs of this bin present in the pathway </li>                
                <li>The column "All orthologs" is the list of the orthologs of all bins present in the pathway </li>
            </ul></li>
            <li> Figure detail
            <ul>
                <li>The Figures in the links: <ul> 
                    <li>The orthologs present from the bin are in <font color="#e7bcd4">▉▉</font></li>
            <li>Troubleshooting
            <ul>
                <li>When the link of the pathway is not loading or not showing anything, it means that there is too much orthologs to show on the figure.
            Try to strip everything after "https://www.kegg.jp/kegg-bin/show_pathway?PATWAY_ENTRY_NUMBER/" to still see the pathway</li>
            </ul></li>
        </ol>
        </p></div>
        </div>
        """
        )

        outfile.write("""
        <h1>Pathways</h1>
            
            <div class="tg-wrap">
                <table class="tg">
                    <tr>
                        <th class="header">Pathways</th>
                        <th class="header">Orthologs</th>
                        <th class="header">All orthologs</th>
                    </tr>
                    <tr>
                        <li><font color="#e7bcd4">▉▉</font>Represent the orthologs of the Bin</th>
                        <th class="header">List of the orthologs of the bin</th>
                        <th class="header">List of the orthologs of all bins</th>

                    </tr>
        """)
        for pathway in dict_global_bin[bin_html]:
            request = 'http://rest.kegg.jp/get/'+pathway
            print(request)
            try:
                req = urlopen(request)
                data = req.read().decode()
                data = data.split('\n')[1]
                pathway_name = re.split("NAME\s+", data)[1]
                set_gene = set()
                if dict_global_bin[bin_html][pathway][2] != "":
                    for gene in dict_global_bin[bin_html][pathway][2]:
                        set_gene.add(gene)
                    set_html = set()
                    for gene in set_gene:
                        set_html.add(gene)
                    list_html = "".join(set_html)
                outfile.write(f"""
                    <tr>
                    <td class="pathway_gene"><a href="https://www.kegg.jp/kegg-bin/show_pathway?{pathway}/{list_html}/default%3d%23e7bcd4">{pathway_name}</a></td>
                    <td class="listgene">
                    """)
                for gene in set_gene:
                    request_gene = 'http://rest.kegg.jp/get/'+gene
                    print(request_gene)
                    try:
                        req_gene = urlopen(request_gene)
                        data_gene = req_gene.read().decode()
                        data_gene = data_gene.split('\n')[1]
                        gene_name = re.split("NAME\s+", data_gene)[1]
                        outfile.write(
                            f"""<a href="https://www.kegg.jp/dbget-bin/www_bget?{gene}">[{gene_name}]</a>; 
                            """)
                    except:
                        outfile.write(f"""[{gene} unknown by KEGG DB]""")

                outfile.write("""</td>
                    <td class="listgene">
                    """)
                for gene in dictgeneral[pathway]:
                    request_gene = 'http://rest.kegg.jp/get/'+gene
                    print(request_gene)
                    try:
                        req_gene = urlopen(request_gene)
                        data_gene = req_gene.read().decode()
                        data_gene = data_gene.split('\n')[1]
                        gene_name = re.split("NAME\s+", data_gene)[1]
                        outfile.write(
                            f"""<a href="https://www.kegg.jp/dbget-bin/www_bget?{gene}">[{gene_name}]</a>; """)
                    except:
                        outfile.write(f"""[{gene} unknown by KEGG DB]""")
                outfile.write("""</td>
                    </tr>
                    """
                    )
            except:
                outfile.write(f"""<tr>
                    <td class="pathway_gene"> {pathway} unknow by the KEGG DATABASE </td></tr> """)
        outfile.close()


def parse(args):
    """
    Main function of the parser

    Arguments: args (object): all the dictionnary argument from argparse
    """
    logger = logging.getLogger(__name__)
    output = args.output
    bins = args.bins
    # try:
    os.makedirs(args.output)

    globalpathwaylist = []
    binnamelist = []
    dict_global_sample = defaultdict(dict)
    dict_global_bin = defaultdict(dict)
    dictgeneral={}
    dictgeneral, dict_global_sample, dict_global_bin, globalpathwaylist, binnamelist = bin_parse(bins,
                                                                                    globalpathwaylist, binnamelist, dictgeneral,
                                                                                    dict_global_sample, dict_global_bin)
    write_html_sample(dictgeneral, dict_global_sample, output,
                      globalpathwaylist, binnamelist)
    write_html_bins(dictgeneral, dict_global_bin, output,
                    globalpathwaylist)

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
