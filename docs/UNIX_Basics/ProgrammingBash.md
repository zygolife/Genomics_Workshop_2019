# Shell programming in BASH

The shell is powerful environment and can be used to encode lots of logic.  At some point it can be useful to switch to even more expressive programming and scripting language like

# Logical operators

**if / elif / else**

```bash
if [ ! -f file.txt ]; then
    echo "Gene X" > file.txt
fi
```

File test operators
 * `-f` - does file exist
 * `-s` - does file exist and is not empty
 * `-d` - does directory exist

In Bash can use also nice operators with the double `[[ ]]`

```bash
if [[ -d testdir && -f testdir/genome.fa ]]; then
  echo "genome is downloaded"
else
    mkdir testdir
    curl -o testdir/genome.fa https://fungidb.org/common/downloads/Current_Release/CimmitisRS/fasta/data/FungiDB-45_CimmitisRS_Genome.fasta
fi
```

Can also set up multiple conditionals

```bash

if [[ $target == "speciesA" ]]; then
  grep $target inputfile > speciesA_hits.txt
elif [[ $target == "speciesB" ]]; then
  grep $target inputfile > speciesB_hits.txt
else
  echo "unknown"
fi
```

# Iterators

Often we want to do something across all files. But the parameters will depend on the filename for input and output. Here are some simple examples.

Iterate through all files ending in ".fa" and run BLAST writing out the result file in a file name which is the same as input but the ".fa" is removed.

# For loops

use `for` to operate on a set of items.

```bash
for file in *.fasta
do
  wc -l $file > $file.lines
done
```

```bash
module load ncbi-blast
for file in $(ls *.fa)
do
  base=$(basename $file .fa) # strip the .fa from the end of the filename
  blastp -query $file -db swissprot -out $base.BLASTP.tab -outfmt 6
done
```

Where we need to provide a number in the teration the command `seq` will generate a sequence of numbers.

```bash
$ seq 5
1
2
3
4
5

$ seq 15 18
15
16
17
18

$ seq 1 2 10 # odd numbers, count by 2's
1
3
5
7
9
$ seq 10 7 # reverse order
10
9
8
7
```

```bash
for n in $(seq 1 10)
do
  bash job.sh $n
done
```

## UNIX Parallel

```bash
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

```bash
module load samtools
parallel -j 10 samtools flagstat {} \> {.}.stats ::: *.bam
```

#
