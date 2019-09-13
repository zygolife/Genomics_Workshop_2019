#!/usr/bin/bash
module load bedtools
# This script will calculate the SNPs per gene by looking at intersection of the two files
ln -s ../../data/genomic/C_lusitaniae.genotypes_all.selected.SNPONLY.lungonly.vcf.gz
ln -s ../../data/genomic/candida_lusitaniae_1.sorted.gff3
# make a gene - only feature file
grep -P "\tgene\t" candida_lusitaniae_1.sorted.gff3 > candida_lusitaniae_1.genes.gff3

bedtools intersect -a candida_lusitaniae_1.genes.gff3 -b C_lusitaniae.genotypes_all.selected.SNPONLY.lungonly.vcf.gz -wb | cut -f9 | perl -p -e 's/ID.+Name=//' | sort | uniq -c | sort -rn > genes_by_snpcount.txt

