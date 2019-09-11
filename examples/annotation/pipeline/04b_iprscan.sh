#!/bin/bash
#SBATCH --ntasks 32 --nodes 1 --mem 96G -p intel 
#SBATCH --time 48:00:00 --out logs/iprscan.%a.%A.log

module unload miniconda2
module load miniconda3
module load funannotate/git-live
module load iprscan
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=strains.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
	if [ ! -d $OUTDIR/${BASE} ]; then
		echo "No annotation dir for ${BASE} did you run 01_predict.sh $N?"
		exit
 	fi
	mkdir -p $OUTDIR/${BASE}/annotate_misc
	XML=$OUTDIR/${BASE}/annotate_misc/iprscan.xml
	IPRPATH=$(which interproscan.sh)
	if [ ! -f $XML ]; then
	    funannotate iprscan -i $OUTDIR/${BASE} -o $XML -m local -c $CPU --iprscan_path $IPRPATH
	fi
done
