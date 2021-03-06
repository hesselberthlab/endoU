#! /usr/bin/env python

import sys
import Bio

#A program to report the 3'dinucleotide for each 5' position in a depth
#file for the MHV antigenomic RNA

#specify fasta file on command line
fasta = file(sys.argv[1])

#read in and parse fasta file, transcribe and reverse complement, then
#reverse again to get the complement

from Bio import SeqIO
for record in SeqIO.parse(fasta, "fasta"):
        #rev_seq = seq[::-1]
        seq = record.seq.reverse_complement().transcribe()
        #seq = record.seq.transcribe()
        rev_seq = seq[::-1]
        #print(rev_seq)
        #seq = seq[::-1]

#specify bedgraph file on the command line
bedgraph = file(sys.argv[2])

#open and parse bedgraph file by line, remove whitespace and split fields
#into strings, name the fields, use the start position to slice
#sequence for the dinucleotide 

for line in (bedgraph):
    if line.startswith('#'): continue 
    fields = [field.strip() for field in line.split('\t')]
    chrom, start, count = fields[:3]

    #fasta sequence starts at 0-coord and depth files start at 1-coord
    position = int(start) - 1

    #slice 1 bp upstream and downstream of start
    dinuc = rev_seq[position-1:position+1]    
    print chrom, '\t', start, '\t', count, '\t', dinuc

