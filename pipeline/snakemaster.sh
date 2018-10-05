#!/usr/bin/env bash

#BSUB -J endoU
#BSUB -o snakemake_%J.out
#BSUB -e snakemake_%J.err
#BSUB -R "select[mem>4] rusage[mem=4] " 

set -o nounset -o pipefail -o errexit -x


args=' 
  -o {log}.out 
  -e {log}.err 
  -J {params.job_name} 
  -R "{params.memory} span[hosts=1] " 
  -n {threads} ' 


snakemake --drmaa "$args" \
    --snakefile Snakefile \
    --jobs 40 \
    --resources all_threads=50 \
    --latency-wait 50 \
    --rerun-incomplete  \
    --configfile config.yaml \
    #--unlock
    #--directory /beevol/home/ancarr/projects/endoU/pipeline
