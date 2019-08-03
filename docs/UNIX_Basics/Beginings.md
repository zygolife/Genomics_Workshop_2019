## Getting started

To login to the cluster we need to use ssh client. This allows secure communication with the cluster. The UCR cluster is accessed using the host `cluster.hpcc.ucr.edu`

```
$ ssh USERNAME@cluster.hpcc.ucr.edu
```

This will start up a UNIX session running on the cluster 'head node'. There are multiple machines which serve as this login node where we can stage our analysis to run on the worker nodes that are on the cluster. Much more detail on the setup of the cluster and resources available at http://hpcc.ucr.edu.