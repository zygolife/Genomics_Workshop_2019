#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 8 --mem 96gb --out logs/metahit.%a.log
# needs to run on intel or short queue Only
module load megahit
MEM=64gb
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
IN=input
OUTDIR=asm_metahit
mkdir -p $OUT
SAMPLEFILE=samples.csv
IFS=,
sed -n ${N}p $SAMPLEFILE | while read SAMPLE NAME
do
  LEFT=$IN/${SAMPLE}_1.fastq.gz
  RIGHT=$IN/${SAMPLE}_2.fastq.gz
megahit --tmp-dir /scratch -1 $LEFT -2 $RIGHT \
-o $OUTDIR/$NAME --mem-flag 1  --memory 0.5 -t $CPU

# or can run with different options like extra sensitivity
#megahit --tmp-dir /scratch -1 left.fq.gz -2 right.fq.gz  --presets meta-
#sensitive -o megahit_sample_sensitive
done
