#!/usr/bin/env bash

#BSUB -J TRS[1-46]
#BSUB -o TRS.%J.%I.out
#BSUB -e TRS.%J.%I.err
#BSUB -n 20
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x


#program to run python script for parsing TRS intervals from bam files

bam=/beevol/home/ancarr/projects/endoU_working/2019-05-29/pipeline/results/bams
TRS_reads=/beevol/home/ancarr/projects/endoU_working/2019-05-29/TRS_analysis/TRS_counts

#b=B6_mock_12_rep1


#list of bam files 

bam_names=(  )

b=${bam_names[$(( $LSB_JOBINDEX -1 ))]}

python TRS_counting.py $bam/${b}_Aligned.sortedByCoord.out.bam > $TRS_reads/$b.txt














    





