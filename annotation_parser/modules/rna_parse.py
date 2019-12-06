#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv


def rna_parse(rna, dict_transcript, dictrna, rna_pathway_list):
        pathway_list=[]
        with open(rna) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                line_count=0
                for row in csv_reader:
                    if line_count==0:
                        line_count+=1
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
                            elif pathway=="":
                                pass
                            else:
                                linenumber = line_count+1
                                print(
                                    f"error unknown pathway {pathway} in line number {linenumber} in {rna}")
                        line_count+=1
        rna_pathway_list = list(set(pathway_list))
        return dict_transcript, dictrna, rna_pathway_list
