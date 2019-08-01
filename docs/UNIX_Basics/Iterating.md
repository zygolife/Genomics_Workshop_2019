## Loops

Often we want to apply the same step to multiple files. Here are some
simple examples.


## For loops

```
for file in *.fasta
do
  wc -l $file > $file.lines
done
```

## UNIX Parallel

```
module load parallel # on HPCC
# sometimes you need to deal with some perl issues on UCRHPCC
# uncomment this next line if problems
# module unload perl

parallel -j 10 grep -c ">" {} \> {}.reads  ::: *.fasta
```

Or more practically if you want to get alignment stats on
many BAM files. This will run alignment stats (`samtools flagstat`) on
a series of files and create files with the .bam part removed.
The parallel special replacement code `{.}` means strip off everything
after the last '.' in the filename.

```
module load samtools
parallel -j 10 samtools flagstat {} > {.}.stats ::: *.bam
```

### Iterating where we need to provide a number

The command `seq` will generate a sequence of numbers

```bash
$ seq 10
1
2
3
4
5
6
7
8
9
10

$ seq 15 20
15
16
17
18
19
20
$ seq 1 2 10 # odd numbers
1
3
5
7
9
$ seq 10 7
10
9
8
7
```
