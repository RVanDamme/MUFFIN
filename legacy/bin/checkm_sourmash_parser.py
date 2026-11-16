#!/usr/bin/env python

import os
import re
import csv
import sys
import shutil
import logging
import argparse


def checkm_parser(dict_checkm_sourmash,checkm):
    with open(checkm) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            # line=row[0].split(",")
            dict_checkm_sourmash[row[0]] = list(row[0:])
        # print(dict_checkm_sourmash)
        return dict_checkm_sourmash

def sourmash_parser(dict_checkm_sourmash,sourmash):
    with open(sourmash) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            dict_checkm_sourmash[row[0]].extend( row[1:] )
        # print(dict_checkm_sourmash)
        return dict_checkm_sourmash


def out_writing(dict_checkm_sourmash):
    with open("classify_step_summary.csv", mode='w', newline='') as out_file:
        out_writer = csv.writer(out_file, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)

        out_writer.writerow(
            ['Bin Id', 'Marker lineage', 'UID', 'genomes', 'markers', 'marker sets',
             '0', '1', '2', '3', '4', '5+', 'Completeness', 'Contamination',
              'Strain heterogeneity', 'GTDB status', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
              )
        for bins in dict_checkm_sourmash:
            out_writer.writerow(dict_checkm_sourmash[bins])
def parse(args):
    """
    Main function of the parser

    Arguments: args (object): all the dictionnary argument from argparse
    """
    logger = logging.getLogger(__name__)
    checkm = args.checkm
    sourmash = args.sourmash

    dict_checkm_sourmash = {}
    dict_checkm_sourmash = checkm_parser(dict_checkm_sourmash,checkm)
    dict_checkm_sourmash = sourmash_parser(dict_checkm_sourmash, sourmash)
    out_writing(dict_checkm_sourmash)


def main():
    # argument handling
    parser = argparse.ArgumentParser(
        description='Parse checkm and sourmash files to output a nice csv summary file')
    parser_input = parser.add_argument_group("required input")

    parser_input.add_argument(
                "-c",
                "--checkm",
                required=True,
                dest="checkm",
                help="checkm stripped file (see MAFIN module called checkm_sourmash_parser)",
                metavar='checkm.tsv'
                )
    parser_input.add_argument(
                "-s",
                "--sourmash",
                required=True,
                dest="sourmash",
                help="sourmash stripped file (see MAFIN module called checkm_sourmash_parser)",
                metavar='sourmash.csv'
                )

    args = parser.parse_args()

    try:
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger("checkm_sourmash_parser")
        parse(args)
    except AttributeError as e:
        logger.debug(e)
        parser.print_help()
        raise


if __name__ == "__main__":
    main()
