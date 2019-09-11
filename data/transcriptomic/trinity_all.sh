#!/usr/bin/bash
#SBATCH -p intel --mem 64gb -N 1 -n 16 --out trinity_all.log
module load trinity-rnaseq
Trinity --left rs_con1_1.fq.gz --right rs_con1_2.fq.gz --CPU 16 --max_memory 64G --seqType fq --output trinity_out_allreads
