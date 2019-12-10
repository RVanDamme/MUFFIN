#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

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
                line_count+=1
        return dictquant
