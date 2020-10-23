#! /usr/bin/env bash

#BSUB -J dinucs[1-2]
#BSUB -o logs/%J.%I.out
#BSUB -e logs/%J.%I.err

#script to run python program to fetch dinucleotides from transcripts
#aligned to the mouse transcriptome and MHV genome 
#in both the -2:-1 and -1:+1 registers 

set -o nounset -o pipefail -o errexit -x

refs=(mrna virus)
fastas=(gencode.vM22.transcripts.fa.gz mhv.mod.fasta.gz)

reps=(rep1 rep2)
offset_types=(minus-one-plus-one zero-plus-two)
offset_coords=("-1|1" "0|2")

ref=${refs[$(($LSB_JOBINDEX - 1))]}
fasta=${fastas[$(($LSB_JOBINDEX - 1))]}

for i in ${!offset_types[@]}; do

    offset_type=${offset_types[$i]}

    offset_coord=${offset_coords[$i]}
    up=$(echo $offset_coord | cut -f1 -d"|")
    down=$(echo $offset_coord | cut -f2 -d"|")

    for rep in ${reps[@]}; do
        results=results/$ref/$rep/$offset_type
        mkdir -p $results
        for tab in $(ls -d ref/$ref/$rep/*); do
            outfile="$results/$(basename $tab .bg.gz).tab.gz"
            python3 fetch_dinucs.py $fasta $tab --up $up --down $down \
                | gzip -c > $outfile
        done
    done
done

