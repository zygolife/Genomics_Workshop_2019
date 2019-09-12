#!/usr/bin/bash
#SBATCH --ntasks 18 --mem 16G --time 2:00:00 -p short

module load hmmer/3
module load python/3
if [ ! -f config.txt ]; then
	echo "Need config.txt for PHYling"
	exit
fi

source config.txt
if [ ! -z $PREFIX ]; then
	rm -rf aln/$PREFIX
fi
# probably should check to see if allseq is newer than newest file in the folder?

./PHYling_unified/PHYling init
./PHYling_unified/PHYling search
./PHYling_unified/PHYling aln -c
pushd phylo
sbatch --time 24:00:00 -p batch fast_run.sh
