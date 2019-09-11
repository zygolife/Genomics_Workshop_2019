#!/usr/bin/bash
#SBATCH -p short --mem 16gb --out logs/DESeq.log
mkdir -p plots reports
Rscript Rscripts/kallisto_DESeq.R > logs/DESeq_kallisto.log
