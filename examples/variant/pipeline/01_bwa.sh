#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out logs/bwa.%a.log

CPU=8

module load bwa
module load samtools

INPUT=input
OUTPUT=aln
SAMPLEFILE=strains.csv

mkdir -p $OUTPUT

IFS=,
sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN
do
	if [ ! -f $OUTPUT/$STRAIN.sam ]; then
		bwa mem -t $CPU genome/Scerevisiae.fasta \
		$INPUT/${BASE}_1.fastq.gz $INPUT/${BASE}_2.fastq.gz > $OUTPUT/$STRAIN.sam 
	fi
	if [ ! -f $OUTPUT/$STRAIN.bam ]; then
		samtools sort -O bam -o $OUTPUT/$STRAIN.bam --threads $CPU -T /scratch $OUTPUT/$STRAIN.sam
	fi
	if [ ! -f $OUTPUT/$STRAIN.bam.bai ]; then
		samtools index $OUTPUT/$STRAIN.bam
	fi
done
