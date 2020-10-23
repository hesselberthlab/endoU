#!/usr/bin/env bash

#BSUB -J vcf[1-46]
#BSUB -o vcf.%J.%I.out
#BSUB -e vcf.%J.%I.err
#BSUB -n 20
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

#Program to filter vcf files for snp positions and quality data 

vcf=/beevol/home/ancarr/projects/endoU_working/2019-05-29/snp_analysis/vcf
fil_vcf=/beevol/home/ancarr/projects/endoU_working/2019-05-29/snp_analysis/vcf/filtered_vcf


#list of vcf files to be filtered

vcf_names=(
RNaseL_ns2_9_rep2
B6_nsp15_9_rep2) 


v=${vcf_names[$(( $LSB_JOBINDEX -1 ))]}


python filter_vcf.py $vcf/$v.var.flt.vcf > $fil_vcf/$v.var.txt







    







