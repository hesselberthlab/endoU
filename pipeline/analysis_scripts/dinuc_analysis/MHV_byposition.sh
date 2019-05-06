#!/usr/bin/env bash

#BSUB -J dinuc[1]
#BSUB -o dinuc.%J.%I.out
#BSUB -e dinuc.%J.%I.err
#BSUB -n 8
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x


bedgraph=B6_MHVV_12.mhv.bg

#b=${bedgraph[$(( $LSB_JOBINDEX -1 ))]}

fasta=/beevol/home/ancarr/data-sets/genome/virus/mhv.mod.fasta

python antisense_byposition_+2.py $fasta $bedgraph | cut -f2,4 > mhv.dinuc.antisense_+2.txt

#python sense_byposition.py $fasta $bedgraph | cut -f2,4 > mhv.dinuc.sense.txt




