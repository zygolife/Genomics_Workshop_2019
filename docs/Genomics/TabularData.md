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


# Tools for working with Tabular data

**BEDTools** The [BEDTools](https://bedtools.readthedocs.io/en/latest/) provides lots of examples

**
