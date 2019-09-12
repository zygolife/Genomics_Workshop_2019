#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out logs/bwa.%a.log --time 2:00:00

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

module load bwa
module load samtools

INPUT=input
OUTPUT=aln
SAMPLEFILE=strains.csv

mkdir -p $OUTPUT

IFS=,
sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN
do
    if [ ! -f $OUTPUT/$STRAIN.bam ]; then
	
	if [ ! -f $OUTPUT/$STRAIN.sam ]; then
	    bwa mem -t $CPU genome/Scerevisiae.fasta \
		$INPUT/${BASE}_1.fastq.gz $INPUT/${BASE}_2.fastq.gz > $OUTPUT/$STRAIN.sam 
	fi
	samtools sort -O bam -o $OUTPUT/$STRAIN.bam --threads $CPU -T /scratch $OUTPUT/$STRAIN.sam
	
	if [ ! -f $OUTPUT/$STRAIN.bam.bai ]; then
	    samtools index $OUTPUT/$STRAIN.bam
	fi
	rm $OUTPUT/$STRAIN.sam
    fi
done
