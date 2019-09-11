#!/bin/bash
#SBATCH --nodes 1 --ntasks 1 -p short --mem 4gb --out logs/init.log

mkdir -p logs
module load kallisto
DB=Rdel_transcripts.fasta
OUT=$(basename $DB .fasta)
for size in 31;   # other values might be 27 21;
do
    if [ ! -f ${OUT}.kallisto_k${size}.idx ]; then
	kallisto index -k $size -i ${OUT}.kallisto_k${size}.idx ${OUT}.fasta
    fi
done

