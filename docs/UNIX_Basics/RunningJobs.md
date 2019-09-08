# Running jobs on the cluster

The HPCC system provides a [useful explanation and tutorial](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html) on running jobs on the system. A job is how to run a set of commands and the machine it is sent to on the cluster will be assigned by the queueing system.

![https://implement.pt/img/hpc_job_sched.png]

# Interactive Job session

Startup an interactive session with 1 cpu, 2gb of RAM. Run on the 'short' queue which gives a maximum of 2 hrs.

```
$ srun -N 1 -n 1 --mem 2gb -p short --pty bash -l
$ module load ncbi-blast
$ curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/swissprot.gz
$ pigz -d swissprot.gz
$ makeblastdb -in swissprot -dbtype prot
```

# Submitting a job

Using an editor (emacs,nano,vi) - create a script that has this content. Call it `makeblastdb.sh`
```
#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 --mem 2gb -out make_swissprot_blast.log
module load ncbi-blast/2.9.0+
curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/swissprot.gz
pigz -d swissprot.gz
makeblastdb -in swissprot -dbtype prot
```

Now submit the job via slurm using the command line
```
$ sbatch makeblastdb.db
```

# Array jobs

A great deal of informatics and genomics can be split into smaller pieces and running independently. This allows taking advantage of parallel job running on high performance computing.

Array jobs is a way to run the same command across a set or scripts where the only thing different is the value of a single variable. This variable is an integer which can be used to determine a particular slice or partition of the data to process.

The simplest way to set these up is to think of array job as a line number in a file. The file is a list of names of other files to process. So for each job run it would process a different line in the file.

This uses the command `--array=1` or to specify a range with `--array=1-10`. These need not be continuous, so `--array=4,7,12-30` is all valid.


# A basic template

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1 # number of CPUS
#SBATCH --mem=8gb # how much memory
#SBATCH --time=0-00:15:00     # 15 minutes
# %A is the job number, %a is the array job number
#SBATCH --output=my.stdout.%A_%a.out
#SBATCH --job-name="Name_of_Job"

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
date
echo "${SLURM_ARRAY_JOB_ID}[${SLURM_ARRAY_TASK_ID}]"
sleep 60
hostname

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
```

Some other SBATCH options include way to set the memory proportional to the number of CPUs allocated.
```
# SBATCH --mem-per-cpu=1G
```
Specify an email to be sent when a job is finished
```
#SBATCH --mail-user=YOUREMAIL
#SBATCH --mail-type=ALL
```
