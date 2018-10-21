#! /usr/bin/env python

''' split a bedgraph by chromosome '''

import sys
import os
import glob

#BG_NAMES = {"mhv" , "tRNA", "r18s" , "r28s" , "r5.8s" , "r5s" , "U6" ,"mRNA"}

bedgraph = file(sys.argv[1])

#with open ('bedgraph', 'rt') as in_file

for line in bedgraph:
    fields = line.split("\t")
    chrom, start, end, count, normalized_count = fields[:5]
    chrom = fields[0].split("|")[0].strip(">")
    if chrom.startswith({bg_name}):
        print chrom, '\t',start, '\t',end,'\t',count,'\t', normalized_count







