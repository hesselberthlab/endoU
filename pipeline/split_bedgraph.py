#! /usr/bin/env python3

''' split a bedgraph by chromosome '''

import sys
import os

chrom_names = ["mhv" , "tRNA", "r18s" , "r28s" , "r5.8s" , "r5s" , "U6", "mRNA"]
name = os.path.basename(sys.argv[1]).split(".")[0]
out_dir = str(sys.argv[2]) 
#sample = str(sys.argv[2]).split(".")[0]

# set up the output files for writing
outputs = {}
for chrom in chrom_names:
    outpath = os.path.join(out_dir, chrom, name + "." + chrom + ".bg") 
    outputs[chrom] = open(outpath, "w")

for line in open(sys.argv[1]):

    chrom, start, end, count, normalized_count = line.split("\t")
    #print(fields)
    file_chrom = chrom.split("|")[0]
    #chrom, start, end, count, normalized_count = fields[0:5]
    # look up the output file by chrom name
    output = outputs[file_chrom]
    
    output.write(chrom + "\t" +  start + "\t" +  end + "\t" + count + "\t" + normalized_count)
    #print(chrom, start, end,count, normalized_count,
          #sep = "\t", file=output)

