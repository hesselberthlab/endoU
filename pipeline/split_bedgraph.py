#! /usr/bin/env python3

''' split a bedgraph by chromosome '''

import sys

chrom_names = ["mhv" , "tRNA", "r18s" , "r28s" , "r5.8s" , "r5s" , "U6", "mRNA"]

# set up the output files for writing
outputs = {}
for name in chrom_names:
    outputs[name] = file("%s.bg" % name, "w")

for line in file(sys.argv[1]):

    fields = line.split("\t")
    chrom = fields[0].split("|")[0].strip(">")
    # look up the output file by chrom name
    output = outputs[chrom]

    print(chrom, start, end,count, normalized_count,
          sep = "\t", file=output)

