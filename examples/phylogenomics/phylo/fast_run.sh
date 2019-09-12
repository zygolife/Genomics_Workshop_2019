#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 32 --mem 24gb --time 8:00:00 -p intel --out fasttree_run.%A.log
module unload miniconda2
module load miniconda3
module unload perl
module load fasttree/2.1.11
NUM=$(wc -l ../prefix.tab | awk '{print $1}')
source ../config.txt
ALN=../$FINALPREF.${NUM}_taxa.$HMM.aa.fasaln
TREE1=$FINALPREF.${NUM}_taxa.$HMM.ft_lg.tre
TREE2=$FINALPREF.${NUM}_taxa.$HMM.ft_lg_long.tre
if [ ! -s $TREE1 ]; then
	FastTreeMP -lg -gamma < $ALN > $TREE1
	echo "ALN is $ALN"
fi
if [ -s $TREE1 ]; then
    perl ../PHYling_unified/util/rename_tree_nodes.pl $TREE1 ../prefix.tab > $TREE2
fi
