#!/usr/bin/bash 

#SBATCH --nodes=1 -p batch
#SBATCH --ntasks=8
#SBATCH  --mem 8gb 
#SBATCH  --time=36:00:00
#SBATCH --job-name genemark
#SBATCH --output=logs/train_Genemarkhmm.%a.%A.out
module unload miniconda3
module load genemarkESET
module load perl/5.20.2

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=annotate

GENERICPARAM=$(realpath lib/Fusarium_euwallaceae.v1.gmes.mod)
mkdir -p $OUTDIR

SAMPLEFILE=strains.csv
BUSCO=sordariomyceta_odb9
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
IFS=,
sed -n ${N}p $SAMPLEFILE | while read BASE SPECIES
do
    strain=$(echo $BASE | perl -p -e 's/_/ /g')
    if [ ! -f $INDIR/$BASE.masked.fasta ]; then
	echo "No genome for $INDIR/$BASE.masked.fasta yet - run 00_mask.sh $N"
	exit
    fi
   
    OUTDEST=$(realpath $OUTDIR/$BASE/predict_misc)
    if [ ! -d $OUTDEST ]; then
	echo "No annotation dir for $OUTDEST - have you already tried to run predict?"
	exit
    fi
    pushd $OUTDEST
    if [ -f gmhmm.mod ]; then
	echo "gmhmm.mod exists for $OUTDIR, my work is done"
	exit
    fi
    mkdir -p genemark
    pushd genemark
    echo "gmes_petap.pl --core $CPU --min_contig 10000 --fungus --ES --sequence ../genome.softmasked.fa"
    nohup gmes_petap.pl --core $CPU --min_contig 10000 --fungus --ES --sequence ../genome.softmasked.fa >& train.log
    if [ -f output/gmhmm.mod ]; then
	rsync -aL  output/gmhmm.mod ../
	echo "Success!"
	popd
    else
	popd
	echo "Failed to make training file"
	echo "Symlinking existing $GENERICPARAM"
	ln -s $GENERICPARAM gmhmm.mod
    fi
    popd
done
