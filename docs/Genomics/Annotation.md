# Genome Annotation

The description of genome annotation will be focused around a particular tool we wrote and which is geared towards fungal genome annotation. [Funannotate](http://funannotate.readthedocs.io) is series of integrated pipelines and analysis steps which can train gene predictors, integrate RNAseq data, and predict and combine best evidence for accurate gene prediction and deposition to GenBank.

Before you start a real project you need to register a BioProject, BioSample for your genome to prepare to deposit your sequences. This will include raw data that went into the assembly. When you register a project to deposit a genome you will get back an email with LOCUS prefix information from NCBI. We will use that locus prefix in our annotation so that we have a prefix for how genes are named.

See example run in `examples/annotation`
## Setup samples.csv

The key columns that are needed are a **locus prefix**, species, and strain, and typically some way to decide what phylum or BUSCO dataset to use.

For Zygo project we use this data file
```text
ProjID,JGISample,JGIProjName,JGIBarcode,SubPhyla,Species,Strain,Notes
1978,1111687,BCCTH,10804.2.180461.GAGCTCA-TTGAGCT,Mortirellomycotina,Mortierella cystojenkinii 1230,1230,
1978,1111688,BCCTN,10804.2.180461.ATAGCGG-ACCGCTA,Mortirellomycotina,Mortierella horticola AD009,AD009,
```

Locus prefix == JGIProjName columns
Species == Species
Strain  == Strain
Clade   == SubPhyla

For Fusarium project we use this format

```text
KOD_792,Fusarium sp. AF-12
NRRL_13338,Fusarium nelsonii
NRRL_13368,Fusarium longipes-1
NRRL_13371,Fusarium buharicum
NRRL_13374,Fusarium longipes-2
```

## Masking

We need to run repeat Masking#!/bin/bash
```bash
#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 8 --nodes 1 --mem 8G --out logs/mask.%a.log
# This script runs Funannotate mask step
# Because this a project focused on population genomics we are assuming the repeat library
# generated for one R.stolonifer is suitable for all to save time this is used
# This expects to be run as slurm array jobs where the number passed into the array corresponds
# to the line in the samples.info file

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

if [ -z $SLURM_JOB_ID ]; then
    SLURM_JOB_ID=$$
fi

INDIR=genomes
OUTDIR=genomes
SAMPLEFILE=strains.csv
N=${SLURM_ARRAY_TASK_ID}
RM_SPECIES=fungi

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
IFS=,

tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
    IN=$(realpath $INDIR/$BASE.sorted.fasta)
    OUT=$(realpath $OUTDIR/$BASE.masked.fasta)

 if [ ! -f $IN ]; then
     echo "Cannot find $BASE.sorted.fasta in $INDIR - may not have been run yet"
     exit
 fi

 if [ ! -f $OUT ]; then

    module load funannotate/git-live
    module unload rmblastn
    module load ncbi-rmblast/2.6.0
    export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config

    mkdir $BASE.mask.$SLURM_JOB_ID
    pushd $BASE.mask.$SLURM_JOB_ID
    funannotate mask --cpus $CPU -i $IN -o $OUT --repeatmasker_species $RM_SPECIES

    mv funannotate-mask.log ../logs/$BASE.funannotate-mask.log
    popd
    rmdir $BASE.mask.$SLURM_JOB_ID
 else
     echo "Skipping ${BASE} as masked already"
 fi
done
```

## RNA-Seq prediction

If we have RNASeq data then funannotate has a step `train` which will train gene predictors using this information.

## Gene prediction

Prediction step
```bash
#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 8 --nodes 1 --mem 24G --out logs/predict.%a.log


module unload python
module unload perl
module unload miniconda2
module load miniconda3
module load funannotate/git-live
module unload ncbi-blast
module load ncbi-blast/2.2.31+

export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/config)
mkdir -p $TEMP

SEED_SPECIES="anidulans"
BUSCOPATH=/srv/projects/db/BUSCO/v9
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=annotate
SEQCENTER=UCR

mkdir -p $OUTDIR

SAMPLEFILE=strains.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
    if [ ! -f $INDIR/$BASE.masked.fasta ]; then
	echo "No genome for $INDIR/$BASE.masked.fasta yet - run 01_mask.sh $N"
	exit
    fi
    # PEPLIB=$(realpath lib/informant.aa) if you want to specify your own peptide libary for extra evidence
    GENOMEFILE=$(realpath $INDIR/$BASE.masked.fasta)
    OUTDEST=$(realpath $OUTDIR/$BASE)
    mkdir $BASE.predict.$SLURM_JOB_ID
    pushd $BASE.predict.$SLURM_JOB_ID

    funannotate predict --cpus $CPU --keep_no_stops --SeqCenter $SEQCENTER \
	--busco_db $BUSCOPATH/$BUSCO --strain "$STRAIN" \
	-i $GENOMEFILE --name $BASE \
	--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta  \
	--min_training_models 100 \
	-s "$SPECIES"  -o $OUTDEST --busco_seed_species $SEED_SPECIES
    popd

    rmdir $BASE.predict.$SLURM_JOB_ID
done
```


## Functional annotation

Run AntiSMASH for secondary metabolite Prediction

```BASH
#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 96G --out logs/antismash.%a.%A.log -J antismash

module load antismash
module unload perl
source activate antismash
which perl

CENTER=UCR
OUTDIR=annotate
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=strains.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
	 if [ ! -d $OUTDIR/${BASE} ]; then
		echo "No annotation dir for ${BASE} did you run 01_predict.sh $N?"
		exit
 	fi
	mkdir -p $OUTDIR/${BASE}/annotate_misc
	antismash --taxon fungi --outputfolder antismash \
	    --asf --full-hmmer --cassis --clusterblast --smcogs --subclusterblast --knownclusterblast -c $CPU \
	    $OUTDIR/predict_results/*.gbk
done
```

Run Interpro scan for domain Prediction
```BASH
#!/bin/bash
#SBATCH --ntasks 32 --nodes 1 --mem 96G -p intel
#SBATCH --time 48:00:00 --out logs/iprscan.%a.%A.log

module unload miniconda2
module load miniconda3
module load funannotate/git-live
module load iprscan
CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

SAMPFILE=strains.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read BASE SPECIES STRAIN BUSCO LOCUS
do
	if [ ! -d $OUTDIR/${BASE} ]; then
		echo "No annotation dir for ${BASE} did you run 01_predict.sh $N?"
		exit
 	fi
	mkdir -p $OUTDIR/${BASE}/annotate_misc
	XML=$OUTDIR/${BASE}/annotate_misc/iprscan.xml
	IPRPATH=$(which interproscan.sh)
	if [ ! -f $XML ]; then
	    funannotate iprscan -i $OUTDIR/${BASE} -o $XML -m local -c $CPU --iprscan_path $IPRPATH
	fi
done
```

## Other functional Prediction
```BASH
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
```

## SBT file

I have put a default SBT file in the lib folder but you should create your own for either each project or a blanket one for all your projects. You can create the SBT with [this site at NCBI](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/).

## Submitting to GenBank

You can submit the produced .sqn file that is in the `annotate_results` folder to NCBI using the [submit site](https://submit.ncbi.nlm.nih.gov/).
