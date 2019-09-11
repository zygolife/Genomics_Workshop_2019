# Transcriptome Assembly With Trinity

[Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) is a de
novo transcript assembly tool for Illumina data. Typically this is
used to generate an assembly of all RNASeq reads into a consensus set
of transcripts. It is able to separate alternatively spliced
isoforms. It is a useful tool for generating a transcript set that can
be used for genome annotation or in absence of a genome, to

If a genome is available Trinity has a running mode call Genome
Guided, which will first align transcripts to the genome and partition
the assembly into sets of reads which are shown to be linked on a
local chromosome region. This improves the accuracy of the assembled
transcripts to better separate paralogs and can speed up the process.

## Data sets

Try the data in examples/transcriptome. The script `trinity_Rs_scf10.sh` has the following steps to generate a de novo assembly of the transcriptome.

```bash
#SBATCH -p intel -N 1 -n 24 --mem 64gb --out trinity.log
module load trinity-rnaseq
MEM=64G
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
Trinity --left Rs_1.fq.gz --right Rs_2.fq.gz --CPU $CPUS --max_memory $MEM --seqType fq --output trinity_out_scf10
```

To identify putative ORFs you can use [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki).

```BASH
module load transdecoder
TransDecoder.LongOrfs -t trinity_out_scf10/Trinity.fasta
```

## Considerations

It is typical to generate large numbers of reads for transcriptomics
since replicate experiments and multiple conditions are
common. Combining these data can lead to very Large datasets and long
running times. These assemblies will also require substantially more
memory. It is useful to try an smaller set of representative
experiments rather than including all replicates in assemblies. There
are approaches to try to normalize the data and apply a _digital
normalization_ approach which reduces the complexity of datasets. Some
of these [steps are built](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Insilico-Normalization) steps are built into Trinity but it may require some
experimentation to get it correct. Also see http://ivory.idyll.org/blog/trinity-in-silico-normalize.html.
