#! /usr/bin/env python

from collections import defaultdict
from Bio import SeqIO
import sys
import pandas as pd

fasta = file(sys.argv[1]) 

dinucleotide_counts = defaultdict(int)

for record in SeqIO.parse(fasta, "fasta"):
    seq = record.seq.transcribe()

for i in range(len(seq) - 1):
    dinucleotide_counts[seq[i:i + 2]] += 1

df = pd.DataFrame(dinucleotide_counts, index=[0])
df2 = pd.melt(df)
print(df2)




