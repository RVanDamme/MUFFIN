#!/usr/bin/env python

import os
import sys
import shutil
import logging
import argparse
from collections import defaultdict

from annotation_parser.modules.version import __version__
from annotation_parser.modules.bin_parse import bin_parse
from annotation_parser.modules.rna_parse import rna_parse
from annotation_parser.modules.quantify_rna import genelevel
from annotation_parser.modules.quant_parse import quant_parse
from annotation_parser.modules.sample_html import write_html_sample
from annotation_parser.modules.bin_html import write_html_bins



def parse(args):
    """
    Main function of the parser

    Arguments: args (object): all the dictionnary argument from argparse
    """
    logger = logging.getLogger(__name__)
    output = args.output
    bins = args.bins
    level = args.level
    rna = args.rna
    # try:
    os.makedirs(args.output)
        
    dict_transcript = {}
    dictrna = {}
    rna_pathway_list=[]
    dict_transcript, dictrna, rna_pathway_list = rna_parse(
        rna, dict_transcript, dictrna, rna_pathway_list)
        
    dictquant = {}
    dictquant = quant_parse(dictquant, level)
    
    genelevel(output, dict_transcript, dictquant)

    globalpathwaylist=[]
    binnamelist=[]
    dict_global_sample=defaultdict(dict)
    dict_global_bin=defaultdict(dict)  
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
    parser = argparse.ArgumentParser(description='Parse annotations files to output nice HTML results')
    parser_input = parser.add_argument_group("required input")
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print version and exit"
    )
    if not '-v' in sys.argv :
        if not '--version' in sys.argv:
            parser_input.add_argument(
                "-b", 
                "--bins", 
                required=True,
                dest="bins", 
                help= "List of bins annotated by eggnog (required)",
                nargs='+', 
                metavar='bin.?.annotation.tsv'
            )
            parser_input.add_argument(
                "-r", 
                "--rna", 
                required=True,
                dest="rna", 
                help= "Annotation file of the rna transcript (required)",
                metavar='sample.annotations.tsv'
            )
            parser_input.add_argument(
                "-l", 
                "--level", 
                required=True,
                dest="level",
                help= "Quantification of the transcripts expression by trinity/salmon (required)",
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
        logger = logging.getLogger(__name__)
        if args.version:
            logger.info(f"capture version {__version__}")
            sys.exit(0)
        parse(args)
    except AttributeError as e:
        logger.debug(e)
        parser.print_help()
        raise
