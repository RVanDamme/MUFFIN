#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import csv


def bin_parse(bins,
          globalpathwaylist, binnamelist,
          dict_global_sample, dict_global_bin, dictrna):
    pathwaylist=[]
    for binfile in bins:
        bin_name=os.path.splitext(os.path.splitext(os.path.basename(binfile))[0])[0]
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
                                            dictko[pathway]=[gene]
                                        else:
                                            dictko[pathway].append(gene)
                                        if pathway in dictrna:
                                            if gene in dictrna[pathway]:
                                                if not pathway in dictactiv:
                                                    dictactiv[pathway]=[gene]
                                                else:
                                                    dictactiv[pathway].append(gene)
                                for module_number in row[10].split(','):
                                    if module_number != '':
                                        if not pathway in dictmodules:
                                            dictmodules[pathway]=[module_number]
                                        else:
                                            dictmodules[pathway].append(module_number)


                            elif pathway == "":
                                pass
                            else:
                                linenumber = line_count+1
                                print(
                                    f"error unknown pathway {pathway} in line number {linenumber} in {binfile}")
                        line_count+=1
        csv_file.close()
        for pathw in dictko.keys():
            number_ko = len(dictko[pathw])
            try:
                number_activ = len(dictactiv[pathw])
            except KeyError:
                number_activ = 0
            try:
                if number_activ != 0:
                    dict_global_sample[pathw][bin_name]=[number_ko, dictko[pathw], number_activ, dictactiv[pathw], dictmodules[pathw]]
                    dict_global_bin[bin_name][pathw] = [number_ko, dictko[pathw], number_activ, dictactiv[pathw], dictmodules[pathw]]
                else:
                    dict_global_sample[pathw][bin_name]=[number_ko, dictko[pathw], number_activ, "", dictmodules[pathw]]
                    dict_global_bin[bin_name][pathw] = [number_ko, dictko[pathw], number_activ, "", dictmodules[pathw]]
            except KeyError:
                if number_activ!=0:
                    dict_global_sample[pathw][bin_name]=[number_ko, dictko[pathw], number_activ, dictactiv[pathw], ""]
                    dict_global_bin[bin_name][pathw] = [number_ko, dictko[pathw], number_activ, dictactiv[pathw], ""]
                else:
                    dict_global_sample[pathw][bin_name]=[number_ko, dictko[pathw], number_activ, "", ""]
                    dict_global_bin[bin_name][pathw] = [number_ko, dictko[pathw], number_activ, "", ""]
    globalpathwaylist = list(set(pathwaylist))
    return dict_global_sample, dict_global_bin, globalpathwaylist, binnamelist
