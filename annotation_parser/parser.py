#!/usr/bin/env python

import sys
import os
import argparse
import csv


# argument handling
parser = argparse.ArgumentParser(description='Parse annotations files to output nice HTML results')

parser_input = parser.add_argument_group("required input")

parser_input.add_argument("-b", "--bins", required = True, dest="bins", help= "List of bins annotated by eggnog (required)",
                            nargs='+', metavar='bin.?.annotation.tsv')
parser_input.add_argument("-r", "--rna", required = True, dest="rna", help= "Annotation file of the rna transcript (required)",
                            metavar='sample.annotations.tsv')
parser_input.add_argument("-l", "--level", required = True, dest="level",
                            help= "Quantification of the transcripts expression by trinity/salmon (required)",
                            metavar='sample_transcript_quant.sf')

parser_output = parser.add_argument_group("output directory")

parser_output.add_argument("-o", "--output", nargs='?', const=1, default="./parser_results",
                           dest="output", metavar="Out_dir_name",
                           help="Name for the output directory")

args = parser.parse_args()

print(args)