#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16 --mem 16gb
#SBATCH --output=logs/annotfunc.%a.%A.log
#SBATCH --time=2-0:00:00
#SBATCH -p intel -J annotfunc
module load funannotate/git-live
module load phobius
CPUS=$SLURM_CPUS_ON_NODE
SAMPFILE=strains.csv
BUSCOPATH=/srv/projects/db/BUSCO/v9
if [ ! $CPUS ]; then
 CPUS=1
fi
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
TEMPLATE=$(realpath lib/zygo_genomes_template.sbt)
MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi
echo "N is $N"
IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
	if [ ! -d $OUTDIR/${BASE} ]; then
		echo "No annotation dir for ${BASE} did you run 01_predict.sh $N?"
		exit
 	fi
	MOREFEATURE=""
	if [[ ! -z $TEMPLATE ]]; then
		 MOREFEATURE="--sbt $TEMPLATE"
	fi
	CMD="funannotate annotate --busco_db $BUSCOPATH/$BUSCO -i $ODIR --species \"$SPECIES\" --strain \"$STRAIN\" --cpus $CPUS $EXTRAANNOT $MOREFEATURE"
	echo $CMD
	eval $CMD
done
