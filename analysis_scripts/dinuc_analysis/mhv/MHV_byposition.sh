#!/usr/bin/env bash

#BSUB -J dinuc[1]
#BSUB -o dinuc.%J.%I.out
#BSUB -e dinuc.%J.%I.err
#BSUB -n 8
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

#script to run python scripts to identify dinucleotide frquencies in MHV
#sense and antisense cyclic phosphate libraries

#requires bedgraph file for input with position information for every base 
#in the viral genome and a viral fasta file


bedgraph=B6_MHVV_12.mhv.bg

#b=${bedgraph[$(( $LSB_JOBINDEX -1 ))]}

fasta=/beevol/home/ancarr/data-sets/genome/virus/mhv.mod.fasta

python sense_dinuc.py $fasta $bedgraph | cut -f2,4 > mhv.sense.txt

python antisense_dinuc.py $fasta $bedgraph | cut -f2,4 > mhv.antisense.txt




