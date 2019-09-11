#!/usr/bin/bash
#SBATCH -p intel --mem 24gb -N 1 -n 16 --out trinity_scf10.log
module load trinity-rnaseq
Trinity --left input/Rs_scf10_1.fq.gz --right input/Rs_scf10_2.fq.gz --CPU 16 --max_memory 24G --seqType fq --output trinity_out_scf10
