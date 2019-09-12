#!/usr/bin/bash

module load aspera

pushd input
/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
mv SRR*/*.gz .

