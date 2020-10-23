#!/usr/bin/env bash

#BSUB -J vcf[1-46]
#BSUB -o vcf.%J.%I.out
#BSUB -e vcf.%J.%I.err
#BSUB -n 20
#BSUB -R "span[hosts=1]"

set -o nounset -o pipefail -o errexit -x


bam=/beevol/home/ancarr/projects/endoU_working/2019-05-29/pipeline/results/bams
bam_mhv=/beevol/home/ancarr/projects/endoU_working/2019-05-29/snp_analysis/mhv_bams
ref=/beevol/home/ancarr/data-sets/genome/combined_genomes/mhv_mouse/mhv.mod.fasta
fil_bcf=/beevol/home/ancarr/projects/endoU_working/2019-05-29/snp_analysis/vcf

#b=B6_mock_12_rep1

#List of bam files
bam_names=( )


b=${bam_names[$(( $LSB_JOBINDEX -1 ))]}

#Parsing bam files for reads alinging to MHV

samtools view -h $bam/${b}_Aligned.sortedByCoord.out.bam MHVA59 | \
samtools view -bS - > $bam_mhv/$b.mhv.bam

#snp calling using mpileup 

bcftools mpileup -Ou -f $ref $bam_mhv/$b.mhv.bam | \
bcftools call -Ou -mv | \
bcftools filter -s LowQual -e '%QUAL<30 || DP<30' > $fil_bcf/$b.var.flt.vcf












    





