#SBATCH -p short -N 1 -n 2 --mem 4gb
module load aspera
pushd input
/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
ln -s ERR*/*.fastq.gz .
