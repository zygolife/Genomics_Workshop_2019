#!/usr/bin/bash
#SBATCH -p short --mem 32gb -N 1 -n 16 --out logs/snp_call.log

module load bcftools
module load samtools
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

GENOME=genome/Scerevisiae.fasta
ALN=aln
bcftools mpileup --threads $CPU -Ou -f $GENOME $ALN/*.bam  | bcftools call -vmO z -o called_variants.vcf.gz

