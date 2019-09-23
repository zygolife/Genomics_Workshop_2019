# Variant identification

To call SNPs and INDELs or identify variation in copy number of genes or genomic region is an important component of population genetics and genomics.

Definitions
* SNP - Single Nucleotide Polymorphism - identified as differences between a reference genome and an individual usually from genomic sequencing. Though RNASeq and other methods can
* INDEL - INsertion/DELetion
* CNV - Copy Number Variation. Could

# SNP and short INDEL

## setup
Download data into the folder, we are using 3 from yeast 1001 genome project.

```bash
cd examples/variant/input
module load aspera
/bigdata/stajichlab/shared/bin/sra_download.pl -ascp -id $ASPERAKEY sra.txt
ln -s ERR*/*.fastq.gz . # symlink these files into this top folder
```

strain.csv file has the three datasets and corresponding yeast strain ID.
```Text
ERR1308590,Y4
ERR1308591,CEY622
ERR1308592,CBS1489
```
Procedure is to align short reads to reference genome, then capture SNPs with different

# Copy Number Variation

Plotting coverage is one simple way to look for copy number variation.

Process data with mosdepth.
