# Getting started

This workshop will emphasize UNIX skills to support running bioinformatics tools for genomic and evolutionary analysis.


# Using the HPCC Cluster.

For simplicity this tutorial assumes access to the [UCR HPCC](http://hpcc.ucr.edu) system, but otherwise you can run nearly everything on your own unix system that could be Ubuntu, OSX, or CentOS, Debian. The aspects showing how to use the queuing system.

 To login to the cluster we need to use ssh client. This allows secure communication with the cluster. The UCR cluster is accessed using the host `cluster.hpcc.ucr.edu`

```
$ ssh -X USERNAME@cluster.hpcc.ucr.edu
```

This will start up a UNIX session running on the cluster 'head node'. There are multiple machines which serve as this login node where we can stage our analysis to run on the worker nodes that are on the cluster. Much more detail on the setup of the cluster and resources available at http://hpcc.ucr.edu.

The `-X` option tells the system to [forward your X11 connection](https://kb.iu.edu/d/bdnt) which is necessary for running interactive graphics (eg showing an image, running a graphical editor program like emacs)

# The command line interface (CLI)

The command line provides an interface to interact with file system, and to execute commands starting applications. You can monitor processes running on a computer system with tools like `top` and ``

## Directories and files

**ls**  - list the files and folders in a directory. Options include `-l` to list with details (long). `-t` list ordered by time created (time). These can be combined as `ls -ls` or `ls -l -s`. Specify a folder to list other than the current directory with another argument `ls -l data`.  


**mkdir** Create a directory. Give the `-p` option to create an necessary sub folders and also to not give warnings if a folder already exists.

```bash
$ mkdir test
$ mkdir Alpha/Beta/Zeta # will give error
$ mkdir -p Alpha/Beta/Zeta # will not give error
```

**cd** Change directory.

```bash
mkdir dir1
cd dir1 # go into a directory
cd .. # go back a directory
cd /bigdata/stajichlab/shared/projects # go to a folder
cd ZyGoLife # go to another folder
cd - # go back to the last folder
pwd
# will show /bigdata/stajichlab/shared/projects
pushd /srv/projects/db/AAFTF_DB/ # push a folder on the stack
cd ..
pwd
popd # will go back to folder /bigdata/stajichlab/shared/projects
pwd
```

**pwd** Print the current working directory

**rmdir** Remove a folder. Only works if folder is empty

**rm** Remove a file or folder. This is command to be careful with. To delete a folder that contains many folders you can use recursive remove.

***WARNING*** be careful, UNIX does not save backups of files you delete.

```rm -rf folder``` - recursively

## Reading Files

**more** See the contents of a text file, one page at a time. Go to the next page with 'space'. Can search for a specific text with slash (`/`).

```bash
$ more Gene_list.txt
```


**less** See the contents of a text file, one page at a time. Less has *more* options than **more** with arrows which will let you navigate up and down pages and a search option - use the slash (`/`) and then type in a search text it will highlight all the options.

```bash
$ less Gene_list.txt
```


**head** see the first lines in a file. By default this is 10 lines. But you can specify as many as you want with `-n LINES` option.   Useful to get the beginning of a report or see what is the header in a spreadsheet file.
```bash
$ head -n 15 FILE.txt
```

**tail** see the last lines in a file.  By running this command can see the last 10 lines by default. Can specify number of lines with `-n LINES` option.  Useful when looking at a log-file and want to see the last reported messages.
```bash
$ tail -n 12 FILE.txt
```

**cat** - display the contents of a file or files, can use this to concatenate a series of files into one.

```bash
cat file1.txt file2.txt file[789].txt > allfiles.txt
```

**grep** - Find lines in a file that match a pattern.

```bash
grep CM002240 data/genomic/Ncrassa.gff
```

**cp** - Copy file
```bash
cp file1 newfile # copy file
cp -r dir1 dir2  # copy directory recursively
```

**mv** - Rename file
```bash
mv file1.txt file_newname.txt
mv dir2/file1.txt dir3/file_newname.txt
```

**ln** - symlink "ln -s"
```bash
ln -s file1.fasta filename.fa
```

**find** - find files by name or other criteria

Specify starting folder to search and then find by pattern

```BASH
find /srv/project/db/CAZY -name "dbCAN*"
cd /bigdata/stajichlab/shared/projects/ZyGoLife/Bacteria_associates/binned_genomes
find . -name '*.aa.fasta'
```

**man** - Manual pages

```
man find
```

Environment variables.

# Setup your environment with access to the data

Checkout this workshop and data from github. The repository is available at this URL https://github.com/zygolife/Genomics_Workshop_2019 and you can also checkout a copy of this to your local machine or on the cluster.

```bash
$ git clone https://github.com/zygolife/Genomics_Workshop_2019.git
$ cd Genomics_Workshop_2019
# when you want to synchronize
$ git pull
# download some other Datasets
$ cd data
$ bash download_fungidb.sh
```

This will checkout a copy of the data and text for the course. As I update the data or text you can synchronize your copy by issuing `git pull` in the folder that was created.

# Problems

1. Make a series of folders
2. List files in folders
3. Count the number of lines in the file `data/UNIX/Nc3H.expr.tab` is this same as number of lines in `data/UNIX/Nc20H.expr.tab`
4. Practice with pushd / popd
