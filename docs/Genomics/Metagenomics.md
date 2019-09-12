# Metagenomes

Metagenome assemnbly and anlysis of mixed samples of multiple organisms from environment or other sampling.

## Assembly

MetaSPAdes
```bash
#!/usr/bin/bash
#SBATCH -p intel -N 1 --ntasks 24 --mem 100G --out logs/metaspades.%a.log

module load SPAdes

SAMPLES=samples.csv
IN=input
OUTDIR=asm_metaspades
MEM=100
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
IFS=,
sed -n ${N}p $SAMPLES | while read SAMPLE NAME
do
  LEFT=$IN/${SAMPLE}_1.fastq.gz
  RIGHT=$IN/${SAMPLE}_2.fastq.gz
	if [ ! -d $OUTDIR/${SAMPLE} ]; then
		metaspades.py --pe1-1 $LEFT --pe1-2 $RIGHT -o $OUTDIR/${SAMPLE} --threads $CPU --mem $MEM
	else
		echo "Already processed $SAMPLE see $OUTDIR/$SAMPLE"
	fi
done
```
Megahit

```BASH
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
```

## Visualization

### Blob Plots

[BlobPlots](https://blobtools.readme.io/docs/blobplot) are useful ways to show possible contamination or to examine the taxonomy diversity assembly quality.

There are 3 steps to run these.

1. Assign taxonomy identity of each contig. This can be done with NT searches or BLASTX searches. The following uses DIAMOND for a relatively fast screening.
See `pipeline/03_blob_blastx.sh`

```BASH
#SBATCH -p intel --mem 64gb -N 1 -n 32 --out logs/blastx.%a.log

module load diamond

SAMPLES=samples.csv
DB=/srv/projects/db/blobPlotDB/uniprot_ref_proteomes.diamond.dmnd

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
IFS=,
type=megahit
OUT=taxonomy
mkdir -p $OUT
sed -n ${N}p $SAMPLES | while read SAMPLE NAME
do
	ASSEMBLY=asm_${type}/$NAME/final.contigs.fa
	echo $ASSEMBLY
	OUT=$OUT/$NAME.blastx.diamond.tab
	if [[ -f $ASSEMBLY && ! -f $OUT ]]; then
    			diamond blastx --query $ASSEMBLY \
			--db $DB \
			--outfmt 6 --sensitive --max-target-seqs 1 \
			--evalue 1e-25 --threads $CPU --out $OUT
	fi
done
```

2. Next need to generate a BAM file to allow calculation of the coverage of each contig.
See `pipeline/04_blob_make_cov.sh`

```BASH
#!/usr/bin/bash
#SBATCH -N 1 -n 16 -p short --mem 64gb --out logs/make_cov.%a.log

module load bwa
module load samtools/1.9
N=${SLURM_ARRAY_TASK_ID}
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi


if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then
 echo "need to provide a number by --array or cmdline"
 exit
fi
BAMDIR=bam
mkdir -p $BAMDIR

SAMPLES=samples.csv
TEMP=/scratch
IN=input
IFS=,
sed -n ${N}p $SAMPLES | while read SAMPLE NAME
do
	type=megahit
	ASSEMBLY=asm_${type}/$NAME/final.contigs.fa
	BAM=$BAMDIR/$NAME.remap.bam
	READ1=$IN/${SAMPLE}_1.fastq.gz
	READ2=$IN/${SAMPLE}_2.fastq.gz
	if [[ -f $ASSEMBLY && ! -f $BAM ]]; then
		if [ ! -f $ASSEMBLY.bwt ]; then
			bwa index $ASSEMBLY
		fi
	bwa mem -t $CPU $ASSEMBLY $READ1 $READ2 | samtools sort --threads $CPU -T $TEMP -O bam -o $BAM -
	fi
	if [[ -f $BAM && ! -f $BAM.bai ]]; then
		    samtools index $BAM
	fi
done
```

3. These taxonomy and coverage data files are used to construct the Blob Plots.
See `pipeline/05_blob_makeplot.sh`

```bash
#!/usr/bin/bash
#SBATCH -p short --mem 8gb -N 1 -n 1 --out logs/blob.%a.log

module load blobtools/1.1.1
source activate blobtools

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
 N=$1
fi

if [ -z $N ]; then
 echo "need to provide a number by --array or cmdline"
 exit
fi
OUT=blobOut
mkdir -p $OUT

SAMPLES=samples.csv
BAMDIR=bam
TAXDIR=taxonomy
OUT=blobOut
ASM=asm
IFS=,
sed -n ${N}p $SAMPLES | while read SAMPLE READ1 READ2
do
	type=megahit
	mkdir -p $OUT
	ASSEMBLY=${ASM}_$type/${NAME}/final.contigs.fa
	PROTTAX=$TAXDIR/$NAME.blastx.diamond.tab.taxified.out
	BAM=$BAMDIR/${NAME}.remap.bam
	#echo $ASSEMBLY
	#echo $PROTTAX
	#echo $BAM
	if [[ -f $ASSEMBLY && -f $PROTTAX && -f $BAM ]]; then
	    if [ ! -f $OUT/$type/$NAME.AA.blobDB.json ]; then
		blobtools create -i $ASSEMBLY -b $BAM -t $PROTTAX -o $OUT/$NAME.AA
	    fi
	    pushd $OUT
	    if [ -f $NAME.AA.blobDB.json ]; then
		if [ ! -f $NAME.AA.blobDB.table.txt ]; then
		    time blobtools view -r all -i $NAME.AA.blobDB.json
		fi

		for rank in phylum order
		do
		    if [ ! -f $NAME.AA.blobDB.json.bestsum.$rank.p8.span.100.blobplot.read_cov.bam0.png ]; then
			blobtools plot -i $NAME.AA.blobDB.json -r $rank
		    fi
		done
	    fi
	    popd
	fi
done
```

### Classification

We can use fast classification tools like kaiju on raw reads or assembled contigs.

**Kaiju and Krona**

See examples/metagenome/pipeline/02_kaiju.sh
```bash
#!/bin/bash -l

#SBATCH -N 1 -n 10 --mem 100gb -p short
#SBATCH --output=logs/kaiju_out.%a.log
module load kaiju
module load KronaTools

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
		kaiju-addTaxonNames -u -r superkingdom,phylum,class,order,family,genus -t $NODES -n $NAMES -i $OUT/$NAME
.out -o $OUT/$NAME.names.out
	fi
	if [ ! -f $OUT/$NAME.out.krona ]; then
	    kaiju2krona -t $NODES -n $NAMES -i $OUT/$NAME.out -o $OUT/$NAME.out.krona
      ktImportText -o $OUT/$NAME.out.krona.html  $OUT/$NAME.out.krona

	fi
	for clade in phylum family genus
	do
	    if [ ! -f $OUT/$NAME.$clade.summary.tsv ]; then
		kaiju2table -t $NODES -n $NAMES -r $clade -o $OUT/$NAME.$clade.summary.tsv $OUT/$NAME.out
	    fi
	done
 fi
done
```

This will produce a .krona.html file in kaiju_out folder. Open this in a web browser to view Krona plot for the taxonomy.

## Binning
