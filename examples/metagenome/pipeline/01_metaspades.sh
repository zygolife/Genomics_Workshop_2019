#!/usr/bin/bash
#SBATCH -p intel --ntasks 24 --mem 100G --out logs/metaspades.%a.log -N 1

module load SPAdes

SAMPLES=samples.csv
IN=input
OUTDIR=asm_metaspades
MEM=100
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=2
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
IFS=,
sed -n ${N}p $SAMPLES | while read SAMPLE NAME
do
  LEFT=$IN/${SAMPLE}_1.fastq.gz
  RIGHT=$IN/${SAMPLE}_2.fastq.gz
	if [ ! -d $OUTDIR/${NAME} ]; then
		metaspades.py --pe1-1 $LEFT --pe1-2 $RIGHT -o $OUTDIR/${NAME} --threads $CPU --mem $MEM
	else
		echo "Already processed $SAMPLE see $OUTDIR/$SAMPLE"
	fi
done
