## Loops


## For loops

```
for file in *.fasta
do
  wc -l $file > $file.lines
done

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
