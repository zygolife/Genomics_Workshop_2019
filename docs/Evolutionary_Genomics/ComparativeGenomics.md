# Comparative genomics

# Ortholog and Paralog identification


OrthoVenn as a way to view datasets quickly [https://orthovenn2.bioinfotoolkits.net/](https://orthovenn2.bioinfotoolkits.net/).

[OrthoFinder](https://github.com/davidemms/OrthoFinder) is a powerful package for Ortholog identification with a lot of documentation.

# Summarizing Gene Function Content

I have developed a simple Comparative Domain content pipeline available [from github here](https://github.com/stajichlab/Comparative_pipeline).

This will support calculation of Pfam, CAZY (with dbCAN), and MEROPS domains and generate summary tables. There are some associated R code for heatmap generation but this is not included in the package right now.

### Secondary Metabolite Clusters

We run antiSMASH for secondary metabolite prediction, usually from the genbank files produced by funannotate.

```bash

#!/usr/bin/bash
#SBATCH -p intel --time 48:00:00 --out antismash.log -N 1 -n 8

module load antismash
module unload perl
source activate antismash


CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

	antismash --taxon fungi --outputfolder antismash \
	    --asf --full-hmmer --cassis --clusterblast --smcogs --subclusterblast --knownclusterblast -c $CPU \
	    genome_file.gbk
```

