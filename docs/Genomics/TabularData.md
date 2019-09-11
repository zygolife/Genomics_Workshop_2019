# Processing Tabular data

## GFF: Gene Feature Format data.

[GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3) format is tab delimited. There are 9 columns. The 9th column structures the annotation information. There is a specific version - [GFF3](http://gmod.org/wiki/GFF3) which provides a description of how features are related.

There is a nice [graphical representation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) of how format can show relationships of gene and sub-features. This is part of the specification of the GFF3 format with examples.

[GTF](https://genome.ucsc.edu/FAQ/FAQformat.html#format4) format is like GFF but the 9th column is structured

## BED (Browser Extensible Data) format

The [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) was first developed at UC Santa Cruz as part of the UCSC genome browser development efforts.  The format can be as simple as 3 columns which describe a chromosome name, start, and end of a feature.  The data are tab delimited so it is easy to use tools like `cut` and `awk` to process the column-delimited data.

```
chromosome  10  20
```

There are extended version of the format which has additional columns to describe names. The format can be used to describe any data that can be located on linear coordinates.

The four column version has a feature name.
```
chromosome  10  20  GeneX
```

## VCF (Variant Call Format)

This format is used to describe SNPs and INDELs. The [format](http://www.htslib.org/doc/vcf.html) is a text file of tabular data separated by tabs. The columns  describe location, quality, and genotype of a polymorphism in a collection of samples.

# Tools for working with Tabular data

**BEDTools** The [BEDTools](https://bedtools.readthedocs.io/en/latest/) provides lots of examples.

If we wanted to find exons

If we wanted to find the sequence reads that overlap a given region.
The following will report the output from the intersection between the files.

```bash
module load bedtools
bedtools intersect -a genomic/NCU03440.gene.bed -b genomic/ncrassa_exons.bed -wb

CM002237	711893	712513	ID=NCU03440;description=AP-2 complex subunit alpha	CM002237	711893	712513	NCU03440-E4
CM002237	712625	714926	ID=NCU03440;description=AP-2 complex subunit alpha	CM002237	712625	714926	NCU03440-E3
CM002237	714989	715278	ID=NCU03440;description=AP-2 complex subunit alpha	CM002237	714989	715278	NCU03440-E2
CM002237	715364	716221	ID=NCU03440;description=AP-2 complex subunit alpha	CM002237	715364	716221	NCU03440-E1
```

There are many more options and a huge number of options. If you are doing anything to examine genomic ranges, BEDTools often has a tool for the job. It can also examine BAM and CRAM file overlap with BED locations or example.


**Tabix** These tabular formats can be indexed with tabix and compressed with bgzip (both available as part of [htslib](http://www.htslib.org/) package).  To use this to seek out subsets of data you can index the file.
```BASH
module load htslib
bgzip data/genomic/ncrassa_exons.bed
tabix data/genomic/ncrassa_exons.bed.gz
```
Now we can retrieve features that overlap subsets or ranges for example. To get all the exons which fall within this region tabix CM002237 from 200000 to 300000 we can use this:
```bash
module load htslib
tabix data/genomic/ncrassa_exons.bed.gz CM002237:200000-300000
```
will produce
```Text
CM002237	199505	201702	NCU03569-E1
CM002237	228054	229260	NCU03560-E4
CM002237	229322	230044	NCU03560-E3
CM002237	230108	230170	NCU03560-E2
CM002237	230229	230587	NCU03560-E1
CM002237	275717	276561	NCU03549-E4
CM002237	276641	276831	NCU03549-E3
CM002237	276895	277104	NCU03549-E2
CM002237	277201	278127	NCU03549-E1
```
