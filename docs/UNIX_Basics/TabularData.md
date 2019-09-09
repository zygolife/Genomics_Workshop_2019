## Tabular data processing
Data store in text files as tables with simple Tab (TSV) or comma (CSV) delimited format is powerful way to represent and organize data. To access these data we can use several existing tools available in UNIX to work with and access these data.

Useful resources
[Unix Data Tools](https://www.oreilly.com/library/view/bioinformatics-data-skills/9781449367480/ch07.html#chapter-07)

## sort to sort data
[Sort notes](https://biodataprog.github.io/2018_programming-intro/Lectures/03_UNIX_DataProcessing.html#3)
The command `sort` is used to sort data numerically. It can sort by different columns in a multi-column file. [Sort manual](http://man7.org/linux/man-pages/man1/sort.1.html)

To Sort a list by alphanumeric
```bash
sort data/UNIX/Athaliana_chr1_genenames.txt | head
```

```bash
sort -r data/UNIX/Athaliana_chr1_genenames.txt | head
```

```sort -n``` - numeric sort

* -d/--dictionary_order : consider only blanks & alphanumeric characters
* -n/--numeric-sort : compare according to string numerical value
* -f/--ignore-case : upper/lower doesn't matter
* -r/--reverse : reverse the order
* -k specify the key positions to sort by

```bash
sort data/UNIX/numbers_only.dat | head -n 10
10
10
12
25
30
34
39
42
49
49
```

```bash
sort -n data/UNIX/numbers_only.dat | head -n 10
7
7
7
10
10
12
25
30
34
39
```
```bash
sort -r -n numbers_floating.dat  | head -n 10
49.6859213710444
49.6454233452118
49.5141651980655
49.2878027550901
48.5007601226085
45.15231125553
43.0392927946809
41.8950131857132
41.7844270115886
39.63172550297467
```
**cut** retrieve column(s) from delimited dataset.

`cut -f1,3` will return columns 1 and 3

The delimiter is assumed to be tab ("`\t`") but can be any specified single character. For example to separate columns by `,` use:
```bash
$ cut -d- -f1 data/UNIX/all_worm_gene_names.dat
```

**uniq** return a unique list from a sorted list which may contain redundancy. Either a filename is supplied or data is passed in through stdin.

`uniq -c` - take a data list and find only those which are unique, indicate a count of the number of times each value appears.

```bash
$ cut -f1 data/UNIX/Ncrassa_OR74A_InterproDomains.tab | uniq | wc -l
6770

$ cut -f1 data/UNIX/Ncrassa_OR74A_InterproDomains.tab | sort | uniq -c | head -n
   2 NCU00004
   2 NCU00005
   4 NCU00006
   1 NCU00007
   5 NCU00008
   7 NCU00010
  11 NCU00018
   7 NCU00019
```

Use `cut` to grab a column. Then por

```bash
cut -f1 data/UNIX/Ncrassa_OR74A_InterproDomains.tab | uniq -c | sort -nr | head
  63 NCU08377
  50 NCU06272
  44 NCU07119
  39 NCU05316
  38 NCU03545
  37 NCU03116
  34 NCU09438
  34 NCU08933
  34 NCU05545
```

In this case we will use a different delimiter `-` as the gene names are written `srh-97` to indicate the alleles first identified.

```bash
cut -d- -f1 data/UNIX/all_worm_gene_names.dat | sort | uniq -c | sort -n | tail
 283 nhr
 292 srh
 691 let
2829 Bm
3317 Ppa
3624 Cbn
5087 Cjp
6159 Cre
6261 Cbr
15366 21ur
```

# Multicolumn sort

Going back to sort, it is helpful to be able to sort on multiple columns. Imagine you wanted to sort by chromosome number AND start position on the chromosome for data.

```bash
head data/UNIX/random_exons.csv
Chr5,27781790,27781888
Chr11,14656670,14656778
Chr3,14560358,14560608
Chr6,6105107,6105282
Chr1,1147485,1147562
Chr4,8124753,8125072
Chr3,31190249,31190365
Chr9,399276,402077
Chr12,22130532,22130707
Chr2,6789670,6789840

$ sort -t, -k1,1 -k2,2n data/UNIX/random_exons.csv | tr , "\t"  | head -n 5

Chr1	12152	12435
Chr1	1147485	1147562
Chr1	4249358	4249468
Chr1	18658403	18658693
Chr1	43214981	43215253
```

This doesn't quite work because Chr12 comes before Chr2 because the sorting is done based on a string for the first column. To achieve this we would need to split the Chr and number up and sort by the number but this is a little more involved.

# Awk

The **awk** program allows for some very basic programming on tabular data. There are a few different flavors of the program but the [GNU Awk](https://www.gnu.org/software/gawk/manual/gawk.html) (`gawk`) is probably the best one to focus on. You can run it with `gawk` or `awk` (both commands point to the same program on HPCC).  


Awk is different from cut in that it assumes the delimiter is any white space (though that can be changed with the `FS=","` - FS: File Separator.

Print columns 1 and 3 in white space delimited data, order this by the new 2nd column to see genes grouped by the domains they have.
```bash
awk '{print $1,$3}' data/UNIX/Ncrassa_OR74A_InterproDomains.tab | sort -k2
```

# Random numbers

Sometimes we need to generate random numbers. This can be done with different commands. Shuff will shuffle a list and then report a certain number of these, which is an effective way to get a random number betwen 0 and 1000 in the following example.

```bash
shuf -i 0-1000 -n 10
120
280
110
324
621
211
851
70
755
108
```


# Problems

1. What are the top 5 most highly expressed genes, based on the FPKM column, in `data/UNIX/Nc3H.expr.tab` and `data/UNIX/Nc20H.expr.tab`
2. 
