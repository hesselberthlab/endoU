#!/usr/bin/env bash

#BSUB -J dinuc[1]
#BSUB -o dinuc.%J.%I.out
#BSUB -e dinuc.%J.%I.err
#BSUB -n 8
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x


#Program to run python scripts to identify dinucleotide positions from a
#provded bedgraph file containing all position for a speciifc rRNA
#transcript and fasta file for the speciifc rRNA

bedgraph=B6_MHVV_12.r18s.bg

#b=${bedgraph[$(( $LSB_JOBINDEX -1 ))]}

fasta=/beevol/home/ancarr/data-sets/genome/mm10/rRNA/mm18S.fasta

python dinuc_byposition.py $fasta $bedgraph | cut -f2,4 > r18s.dinuc.txt





