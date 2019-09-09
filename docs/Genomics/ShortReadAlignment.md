# Alignment

Sequence alignment to a reference is critical to incorporating 

## Short Read alignment



## SAM/BAM format

This is a column alignment format. SAM is the uncompressed text version of the file. BAM is the binary version and will be compressed. To view a SAM file can just use `less` or `more` and even `grep` to find things, when files are in BAM format will need to use `samtools`

## CRAM format

[CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)) is a Compressed format for alignment data which uses a reference genome to reduce the size of the files by only showing differences to the reference genome. This can achieve 30-60% disk space savings. Many tools can use this format (samtools, mosdepth, bedtools) but it is imperative you have access to the original reference sequence file in FASTA format otherwise data cannot be reextracted.

## Long read alignment

* Minimap2

# Manipulating SAM/BAM files

Sort a BAM file

```
module load samtools
samtools sort -o sorted.bam unsorted.sam
```

As a job - here is a script called `sort.sh`, you would submit it as `sbatch sort.sh`
```
#!/usr/bin/bash
#SBATCH -p short -N 1 -n 8 --mem 8gb --out sort_sam.log
module load samtools
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
samtools sort --threads $CPU -o sorted.bam -T /scratch -


## SAMTools

**Create a BAM file from SAM**


# Summarize Read Depth

[mosdepth](https://github.com/brentp/mosdepth)
