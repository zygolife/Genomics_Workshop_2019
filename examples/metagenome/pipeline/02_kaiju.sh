#!/bin/bash -l

#SBATCH -N 1 -n 10 --mem 100gb -p short
#SBATCH --output=logs/kaiju_out.%a.log
module load kaiju

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
SAMPLES=samples.csv
DBFOLDER=/opt/linux/centos/7.x/x86_64/pkgs/kaiju/share/DB/
DB=$DBFOLDER/nr/kaiju_db_nr.fmi
NODES=$DBFOLDER/nodes.dmp
NAMES=$DBFOLDER/names.dmp
ASM=asm_megahit
OUT=kaiju_out
mkdir -p $OUT
IFS=,
sed -n ${N}p $SAMPLES | while read SAMPLE NAME
do
	QUERY=$ASM/$NAME/final.contigs.fa
echo "QUERY=$QUERY"
 if [ -f $QUERY ]; then
	 if [ ! -f $OUT/$NAME.out ]; then
		kaiju -v -t $NODES -f $DB -z $CPU -a "greedy" -e 5 -s 200 -i $QUERY -o $OUT/$NAME.out
	fi
	if [ ! -f $OUT/$NAME.names.out ]; then
		kaiju-addTaxonNames -u -r superkingdom,phylum,class,order,family,genus -t $NODES -n $NAMES -i $OUT/$NAME.out -o $OUT/$NAME.names.out
	fi
	if [ ! -f $OUT/$NAME.out.krona ]; then
	    kaiju2krona -t $NODES -n $NAMES -i $OUT/$NAME.out -o $OUT/$NAME.out.krona
	fi
	for clade in phylum family genus
	do
	    if [ ! -f $OUT/$NAME.$clade.summary.tsv ]; then
		kaiju2table -t $NODES -n $NAMES -r $clade -o $OUT/$NAME.$clade.summary.tsv $OUT/$NAME.out
	    fi
	done
 fi
done
