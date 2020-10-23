#! /usr/bin/env python

import sys
import re
from pysam import AlignmentFile
import pandas as pd

#program to identify reads containing leader sequence and TRS sequences
#from bamfiles aligned to the MHV genome

bamfile=AlignmentFile(file(sys.argv[1]), "rb")

pos = dict()

for read in bamfile.fetch("MHVA59"):

    leader = "TTTAAATCTAA"

    if re.search(leader, read.seq):

        CIGAR = read.cigartuples
        
        d = dict()

        [d [t [0]].append(t [1]) if t [0] in list(d.keys())
         else d.update({t [0]: [t [1]]}) for t in CIGAR]

        key = read.query_name

        if (0 in d) and (3 in d):
            pos[read.query_name] = [sum([read.reference_start + 1 + d.get(0)[0]] + d.get(3)), read.reference_end] 

        elif (0 in d):
            pos[read.query_name] = [read.reference_start + 1 , read.reference_end]
           # pos[read.query_name] = sum([read.reference_start + 1 + d.get(0)[0]])

df = pd.DataFrame.from_dict(pos, orient = 'index')
df.columns = ['start_pos' , 'end_pos']

#df.columns = ['position']

with pd.option_context('display.max_rows', None):
    print(df)













