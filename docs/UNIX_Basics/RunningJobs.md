# Running jobs on the cluster

The HPCC system provides a [useful explanation and tutorial](http://hpcc.ucr.edu/manuals_linux-cluster_jobs.html) on running jobs on the system. A job is how to run a set of commands and the machine it is sent to on the cluster will be assigned by the queueing system.

![https://implement.pt/img/hpc_job_sched.png]

# Interactive Job session

Startup an interactive session with 1 cpu, 2gb of RAM. Run on the 'short' queue which gives a maximum of 2 hrs.

```
srun -N 1 -n 1 --mem 2gb -p short --pty bash -l
```