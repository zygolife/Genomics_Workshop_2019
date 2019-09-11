#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2

module load bwa

bwa index genome/Scerevisiae.fasta

