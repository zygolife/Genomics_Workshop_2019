#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem 48gb -p intel
#SBATCH --time=3-00:15:00   
#SBATCH --output=logs/train_annot_02.%A.log
#SBATCH --job-name="TrainFun"
module unload python
module unload perl
module unload miniconda2
module load miniconda3
module load funannotate/git-live
PASAHOMEPATH=$(dirname `which Launch_PASA_pipeline.pl`)
TRINITYHOMEPATH=$(dirname `which Trinity`)
export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
CPUS=$SLURM_CPUS_ON_NODE

MEM=48G

if [ ! $CPUS ]; then
 CPUS=2
fi


if [ ! -f config.txt ]; then
 echo "need a config file for parameters"
 exit
fi
ODIR=funannot

source config.txt
if [ ! $MASKED ]; then 
 echo "NEED TO EDIT CONFIG FILE TO SPECIFY THE INPUT GENOME AS VARIABLE: MASKED=GENOMEFILEFFASTA"
 exit
fi

funannotate train -i $MASKED -o funannotate --PASAHOME $PASAHOMEPATH \
 --TRINITYHOME $TRINITYHOMEPATH \
 -s RNASeq_fwd.fq.gz \
   --stranded no --jaccard_clip --species "coccidioides_immitis" --cpus $CPUS --memory $MEM
