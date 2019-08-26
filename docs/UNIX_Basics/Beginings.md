## Getting started

To login to the cluster we need to use ssh client. This allows secure communication with the cluster. The UCR cluster is accessed using the host `cluster.hpcc.ucr.edu`

```
$ ssh USERNAME@cluster.hpcc.ucr.edu
```

This will start up a UNIX session running on the cluster 'head node'. There are multiple machines which serve as this login node where we can stage our analysis to run on the worker nodes that are on the cluster. Much more detail on the setup of the cluster and resources available at http://hpcc.ucr.edu.

### The command line interface (CLI)

The command line or 
### Directories and files

**ls**  - list the files and folders in a directory. Options include `-l` to list with details (long). `-t` list ordered by time created (time). These can be combined as `ls -ls` or `ls -l -s`. Specify a folder to list other than the current directory with another argument `ls -l data`.  


**mkdir** Create a directory. Give the `-p` option to create an necessary sub folders and also to not give warnings if a folder already exists.

```bash
$ mkdir test
$ mkdir Alpha/Beta/Zeta # will give error
$ mkdir -p Alpha/Beta/Zeta # will not give error
```

***rmdir*** Remove a folder. Only works if folder is empty

***rm*** Remove a file or folder. This is command to be careful with. To delete a folder that contains many folders. 

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




