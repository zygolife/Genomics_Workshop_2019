# RNASeq analysis

To Analyze RNASeq data to estimate expression level we need to run
alignment to either transcripts or genome.

This study has multiple timepoints for RNASeq of spore germination in _Rhizopus delemar_.

https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP148808

# Downloading RNASeq data

Go to examples/transcriptome

Run pipeline/00_initialize.sh
```bash
bash pipeline/00_initialize.sh
```
This will download Rdel_transcripts.fasta from [FungiDB](https://fungidb.org/common/downloads/Current_Release/RdelemarRA99-880/fasta/data/). This is the transcripts for Rhizopus delemar. The script will also setup symlink to existing data on cluster or optionally download additional datasets from SRA. If you want to add more steps to download additional raw data from SRA you can follow [the link](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP148808) to the SRA Project and specify additional accession numbers for downloading.

# Aligning RNASeq data with kallisto

First need to index Transcript database. This is accomplished in script
`pipeline/01_index.sh`

```bash
#!/bin/bash
#SBATCH --nodes 1 --ntasks 1 -p short --mem 4gb --out logs/init.log
mkdir -p logs
module load kallisto
DB=Rdel_transcripts.fasta
OUT=$(basename $DB .fasta)
for size in 31;   # other values might be 27 21;
do
    if [ ! -f ${OUT}.kallisto_k${size}.idx ]; then
	kallisto index -k $size -i ${OUT}.kallisto_k${size}.idx ${OUT}.fasta
    fi
done
```

Next the reads can be aligned to these transcripts. We can use array jobs to run this. There are 6 datasets in the samples.csv (first line is a header).

```Text
SampleId,Condition,Replicate,Read1,Read2
24p,24hr,1,SRR7208623_1.fastq.gz,SRR7208623_2.fastq.gz
24a,24hr,2,SRR7208621_1.fastq.gz,SRR7208621_2.fastq.gz
24,24hr,3,SRR7208619_1.fastq.gz,SRR7208619_2.fastq.gz
3p,3hr,1,SRR7208587_1.fastq.gz,SRR7208587_2.fastq.gz
3a,3hr,2,SRR7208584_1.fastq.gz,SRR7208584_2.fastq.gz
3,3hr,3,SRR7208583_1.fastq.gz,SRR7208583_2.fastq.gz
```

We will run arrayjobs to align these transcripts with kallisto and generate read counts per transcript.

```BASH
sbatch --array=1-6 pipeline/02_kallisto.sh
```
This code
```bash
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
```

When all are finished (check with `squeue -u $USER`). Then you can run the R script to generate heatmap and some summary diff expression.

```BASH
sbatch pipeline/03_DESeq.sh
```
This code is
```bash
#!/usr/bin/bash
#SBATCH -p short --mem 16gb --out logs/DESeq.log
mkdir -p plots reports
Rscript Rscripts/kallisto_DESeq.R > logs/DESeq_kallisto.log
```

## Visualizating and interacting with data

[Degust](http://degust.erc.monash.edu/degust/)
