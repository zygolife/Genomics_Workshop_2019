#!/bin/bash
#SBATCH --nodes 1 --ntasks 8 --mem 4G --time 2:00:00 --out logs/Rhizopus_kallisto.%a.log -p short

module load kallisto
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}
INDIR=input
OUTDIR=results/kallisto
SAMPLEFILE=samples.csv
# maybe fix this so all prefixes are consistent!
TRANSCRIPTS=Rdel_transcripts.kallisto_k31.idx

mkdir -p $OUTDIR

if [ -z $N ]; then
 N=$1
 if [ -z $N ]; then
     echo "cannot run without a number provided either cmdline or --array in sbatch"
     exit
 fi
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE CONDITION REP READ1 READ2

do
 OUTNAME=${CONDITION}.r${REP}
 if [ ! -f $OUTDIR/$OUTNAME/abundance.h5 ]; then
     kallisto quant -i $TRANSCRIPTS -o $OUTDIR/$OUTNAME -t $CPU --bias $INDIR/$READ1 $INDIR/$READ2
 fi
done
