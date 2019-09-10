# Genome Assembly

We will perform a genome assembly for an _Aspergillus fumigatus_ strain which I've previously assembly using public data.

We will get this data from a public dataset stored in NCBI. First let's look at the BioProject page [PRJDB3064](http://ncbi.nlm.nih.gov/bioproject/PRJDB3064). This page provides links to the 8 patient isolated cultures.

One strain we will focus on is stored as two SRA datasets [DRR022919](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRR022919) [DRR022920](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DRR022920).

This study has two files for the same strain, one only has 268Mb of data (DRR022919), the other has 1.3Gb (DRR022920). It is probably sufficient to run the assembly process with only the DRR022920 file to keep things simple but we can discuss how you would approach this if you had 2 files.

# Obtaining the sequence data

We are going to download the data using the SRA toolkit and also direct download from EBI of already formatted FASTQ files. The EBI site takes a little longer to download from but provides already formatted FASTQ. The NCBI site is faster to download from usually but you still have to convert the download from the .sra format to .fastq which can take time.  Below I show you several

If you running this NOT on the biocluster you can download it using curl or aspera from the EBI. You can also take advantage of the fact that I have already downloaded the data into the folder `/bigdata/stajichlab/shared/tmp/Afum` so you can simply copy of symlink to this folder.

## Downloading with Curl

`curl`, `wget`, or `lftp` will all allow you run command line download from EBI using a slow but usually accessible method on most computers. By going to the EBI website you can search for the SRA accession number and then bring up this page for the record [DRR022920 at EBI](https://www.ebi.ac.uk/ena/data/view/DRR022920&display=html) and [DRR022919 at EBI](https://www.ebi.ac.uk/ena/data/view/DRR022919&display=html).

Lots more information on downloading at the [EBI read download site](https://www.ebi.ac.uk/ena/browse/read-download).

![EBI page](img/EBI_DRR022920.png)

This site has links to FASTQ files you can use to download. I would probably not do this unless the example further below don't work for you since this will be the slowest download speed. Also for simplicity during the tutorial I have made a copy of all the data so these more like examples when you start to define the data you want to process.

```bash
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR022/DRR022920/DRR022920_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR022/DRR022920/DRR022920_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR022/DRR022919/DRR022919_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/DRR022/DRR022919/DRR022919_2.fastq.gz
```

## Download with aspera

[Aspera](https://downloads.asperasoft.com/connect2/) a very fast download interface. You can access it on the HPCC using
```BASH
module load aspera
ascp -QT -l 1000m -P33001 -i $ASPERAKEY era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/DRR022/DRR022920/DRR022919_1.fastq.gz .
ascp -QT -l 1000m -P33001 -i $ASPERAKEY era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/DRR022/DRR022920/DRR022919_2.fastq.gz .
ascp -QT -l 1000m -P33001 -i $ASPERAKEY era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/DRR022/DRR022920/DRR022920_1.fastq.gz .
ascp -QT -l 1000m -P33001 -i $ASPERAKEY era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/DRR022/DRR022920/DRR022920_2.fastq.gz .
```
This requires knowing ahead of time the path to the folder at EBI to run the aspera download which I determined from the FASTQ FTP download.

However there is a nice Perl script I use which simply requires specifying a file with SRA accessions and will download for you using either FTP or aspera.

You can download it [from here](https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl). It is also available on the biocluster in `/bigdata/stajichlab/shared/bin/sra_download.pl`. To run it make a file with one accession per line, I like to call it `sra.txt`. Here is an example
```
DRR022919
DRR022920
```
To run the download script you can run it as follows on HPCC.

```bash
module load aspera
/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
```

# Setting up analysis

To instead symlink to these data I would suggest focusing only on the DRR022920 accessions.

First make a folder in bigdata to begin your analysis.

```bash
cd ~/bigdata/
mkdir -p assembly/input
cd assembly/input
# create a symlink to the datasets
ln -s /bigdata/stajichlab/shared/tmp/Afum/DRR022920/* .
ln -s /bigdata/stajichlab/shared/tmp/Afum/DRR022919/* .
cd ..
```

Because we want to support a general way of running assemblies across any dataset, we will setup this folder and system in a way that will enable running across multiple isolates.  To do this we need a file that keeps track of a few bits of data about each dataset.  One assumption for the time being is that there is ONE set of paired sequence files for a given isolate you want to assemble. It is possible to recode this so that it can be smarter about taking multiple files OR you can do the really simple things and combine data files before running the analysis steps so that there is only one file per isolate.

To start out we need a data file that we can use in the system to specify information about the dataset. I like to call my file samples.csv but you can vary from that easily as you will see. I called it .csv because we will use commas to separate the data in the file, but it is also possible to do this with tab or even white space as the separator. We will start this off with a samples.csv which has 3 columns, the FATQ READS file basename, the sample/strain name, the expected phyla this is from for some of the optional data cleaning steps. The first is the name of the READ file - in our case that will SRA accession number `DRR022920`. Based on the Biosample record at NCBI or EBI this sample is given the strain designate 'IFM 59356-1'. Because white spaces are often issues I like to change all white spaces to `_` for ease. So we will use the name 'IFM_59356-1'. So our samples file will have for now, one line:

```text
DRR022920,IFM_59356-1,Ascomycota
```

I've created a package I call AAFTF to simplify running all the steps. I'll walk through the parts of the script which relate to each of the pipeline steps. We'll call this job script `01_AAFTF.sh` as I like to number the scripts so I can remember what order they should run in. I typically make a folder called `pipeline` where I keep this script. I also store all the log files that run during an analysis in a folder call `logs` to keep the folder from getting cluttered.

```bash
mkdir -p assembly/pipeline assembly/logs
```
The beginnings of this script are going to be listed here and then I will provide the entire [script here](https://github.com/zygolife/Genomics_Workshop_2019/blob/master/examples/assembly/pipeline/01_AAFTF.sh).

```bash
#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 64gb --out logs/AAFTF.%a.log

MEM=64
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

module load AAFTF

FASTQ=input
SAMPLEFILE=samples.csv
ASM=genomes
WORKDIR=working_AAFTF
mkdir -p $ASM $WORKDIR
if [ -z $CPU ]; then
    CPU=1
fi
IFS=, # set the delimiter to be ,
# read is a step which will process delimited data
# each column will need to match a variable name
# so the number of variables after 'read' needs to be
# same number as number of columns
sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN PHYLUM
do
  # initialize some variables based on the input names
  ASMFILE=$ASM/${STRAIN}.spades.fasta
  VECCLEAN=$ASM/${STRAIN}.vecscreen.fasta
  PURGE=$ASM/${STRAIN}.sourpurge.fasta
  CLEANDUP=$ASM/${STRAIN}.rmdup.fasta
  PILON=$ASM/${STRAIN}.pilon.fasta
  SORTED=$ASM/${STRAIN}.sorted.fasta
  STATS=$ASM/${STRAIN}.sorted.stats.txt
  LEFTIN=$FASTQ/${BASE}_1.fastq.gz
  RIHGHTIN=$FASTQ/${BASE}_2.fastq.gz
  LEFTTRIM=$WORKDIR/${BASE}_1P.fastq.gz
  RIGHTTRIM=$WORKDIR/${BASE}_2P.fastq.gz
  LEFT=$WORKDIR/${BASE}_filtered_1.fastq.gz
  RIGHT=$WORKDIR/${BASE}_filtered_2.fastq.gz
  echo "$BASE $STRAIN"
```
## Data filtering

Generally it is a good idea to do some data filtering to remove contaminating sequence reads like PhiX and vector sequences.

This step will trim low quality bases from reads and remove reads that are too short after this respecting also that the data are paired-end. This uses [bbduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/). You can set a minimum length of trimmed read, the default is 50 but you can set it to be lower if you have very short read data with the option  `--minlength 35`
```bash
AAFTF trim --method bbduk --memory $MEM --left $LEFTIN --right $RIGHTIN -c $CPU -o $WORKDIR/${BASE}
```
This step will filter any reads that match known contaminants using bbduk. It will search for PhiX and vector sequence by default. You can also specify genbank accessions of nucleotide sequences you want to additionally filter with '-a', eg '-a NC_010943.1'. Note that specifying a lot of extra datsets will require a lot more memory so be judicious or consider doing this in a series of looped steps.

```bash
  AAFTF filter -c $CPU --memory $MEM -o $WORKDIR/${BASE} --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk
```

## Assembly


This step runs SPAdes assembler. This was decided as the default assembler that worked well ebough for most of our fungal projects. This step could be replaced with other assemblers, the fastq files are already in the folder and in the `$LEFT` and `$RIGHT` variables so if you have a preference to try out other assemblers for Illumina data this step can be replaced.

```bash
AAFTF assemble -c $CPU --left $LEFT --right $RIGHT --memory $MEM \
	      -o $ASMFILE -w $WORKDIR/spades_${STRAIN}
# this part is just cleanup since we don't need these
# folders if assembly completed, saves space
 if [ -s $ASMFILE ]; then
	    rm -rf $WORKDIR/spades_${STRAIN}/K?? $WORKDIR/spades_${STRAIN}/tmp
	fi

	if [ ! -f $ASMFILE ]; then
	    echo "SPADES must have failed, exiting"
	    exit
	fi
```
## Vector cleaning

Even though we searched for vector sequence already some adaptors or other contaminants are sometime still present. We have adopted a search strategy that is identical to the screen NCBI uses when validating genomes for deposition. If a contig has adaptor sequence internally it will be split by this sequence, the adapator removed and the contig split into 2 or more pieces. Adaptor or vector on the terminal ends of a contig are simply removed. This process runs iteratively until all assembled sequences pass this quality control.

```bash
if [ ! -f $VECCLEAN ]; then
AAFTF vecscreen -i $ASMFILE -c $CPU -o $VECCLEAN
fi
```

## Purge  contamination

This step we have found is probaly not usually very necessary except in cases of heavy contamination (and even then it is probably not aggressive enough). But it adopts a rapid screening approach to classify every contig taxonomically based on hits in a pre-built database as part of [sourmash](https://sourmash.readthedocs.io/en/latest/).  It uses the defined phylum to identify contaminants and flag those contigs and remove them from the generated output file. A report is also generated in the same folder which will indicate the classification for every contig. I think this step could use some refinement it seems to need a better or more comprehensive database. We may try to switch to [Kaiju](https://github.com/bioinformatics-centre/kaiju) in future.

```bash
if [ ! -f $PURGE ]; then
AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $LEFT --right $RIGHT
fi
```

# Remove duplicates and Pilon
There are several additional steps commented out - these represent if you want to try and remove any redundancy in the assembly which can occur usually due to transposon or repeat elements that do not assemble well and generate multiple contigs which are nearly identical. The `rmdup` step will perform that.

The pilon step runs iterations of polishing an assembly.

# Sort the contigs largest to smallest
```bash
if [ ! -f $SORTED ]; then
# AAFTF sort -i $PILON -o $SORTED
AAFTF sort -i $VECCLEAN -o $SORTED
fi
```

# Generate summary statistics about the assembly
```bash
if [ ! -f $STATS ]; then
AAFTF assess -i $SORTED -r $STATS
fi
```
