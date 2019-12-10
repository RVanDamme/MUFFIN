#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

def genelevel(output, dict_transcript, dictquant):
    out = output+"/genelevel.csv"
    with open(out, mode='w', newline='') as gene_file:
        gene_write = csv.writer(gene_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        gene_write.writerow(["transcript id","pathways","Number of reads"])
        for transcript in dict_transcript.keys():
            try:
                gene_write.writerow(
                    [transcript, dict_transcript[transcript], dictquant[transcript]]
                    )
            except KeyError:
                pass
    return
