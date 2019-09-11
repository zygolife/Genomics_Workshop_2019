#!/usr/bin/bash
#SBATCH -p short --mem 32gb -N 1 -n 16 --out logs/snp_call.log

module load bcftools
module load samtools

GENOME=genome/Scerevisiae.fasta
ALN=aln
bcftools mpileup -Ou -f $GENOME $ALN/*.bam  | bcftools call -vmO z -o called_variants.vcf.gz

