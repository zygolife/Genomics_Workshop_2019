#!/usr/bin/bash
#SBATCH -N 1 -n 8 --mem 4gb --out logs/gene_tree_raxml_aa.%a.log -p short

module load RAxML
if [ -f config.txt ]; then
	source config.txt
else
	echo "Need config for outgroup and prefix"
	exit
fi

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

OUTDIR=gene_trees
ALN=$ALN_OUTDIR/$HMM
FILEEXT=aa.trim
raxml=raxmlHPC-PTHREADS-AVX2
mkdir -p $OUTDIR

FILE=$(ls $ALN/*.${FILEEXT} | sed -n ${N}p )
if [ ! $FILE ]; then
 echo "No input file - check $SLURM_ARRAY_TASK_ID or input number, N=$N FOLDER=$ALN ext=$FILEEXT"
 exit
fi
PHY=$FILE
base=$(basename $FILE .$FILEEXT)
OUTPHY=${base}.phy
if [ ! -f $OUTPHY ]; then
    perl -p -e 's/>([^\|]+)\|/>$1 /' $PHY > $OUTDIR/$OUTPHY
fi

pushd $OUTDIR
echo $base
if [ ! -f  RAxML_info.$base"_ML_PROTGAMMA" ]; then
    $raxml -T $CPU -o $OUT -# 100 -x 221 -f a -p 31 -m PROTGAMMAAUTO -s $OUTPHY -n ${base}_ML_PROTGAMMA
fi
