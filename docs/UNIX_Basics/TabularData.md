## Tabular data processing
Data store in text files as tables with simple Tab (TSV) or comma (CSV) delimited format is powerful way to represent and organize data. To access these data we can use several existing tools available in UNIX to work with and access these data.

Useful resources
[Unix Data Tools](https://www.oreilly.com/library/view/bioinformatics-data-skills/9781449367480/ch07.html#chapter-07)

[Sort notes](https://biodataprog.github.io/2018_programming-intro/Lectures/03_UNIX_DataProcessing.html#3)
**sort** - sort data numerically, can sort by different columns in a multi-column file. [manual](http://man7.org/linux/man-pages/man1/sort.1.html)

```sort filename ```

```sort -r``` - reverse ordered list

```sort -n``` - numeric sort

* -d/--dictionary_order : consider only blanks & alphanumeric characters
* -n/--numeric-sort : compare according to string numerical value
* -f/--ignore-case : upper/lower doesn't matter
* -r/--reverse : reverse the order
* -k specify the key positions to sort by

```
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

```
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

**uniq**

```uniq -c``` - take a data list and find only

**cut**

```cut -f1,3```

The delimiter is assumed to be tab ("`\t`") but can be any specified single character. For example to separate columns by `,` use:
```
cut -d, -f1,3
```

### Awk

The **awk** program allows for some very basic programming on tabular data. There are a few different flavors of the program but the [GNU Awk](https://www.gnu.org/software/gawk/manual/gawk.html) (`gawk`) is probably the best one to focus on. You can run it with `gawk` or `awk` (both commands point to the same program on HPCC).  


Awk is different from cut in that it assumes the delimiter is any white space (though that can be changed with the `FS=","` - FS: File Separator.

Print columns 3 and 7 in white space delimited data
`awk '{print $3,$7}' data/UNIX/`


**
