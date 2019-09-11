#!/usr/bin/bash

if [ ! -f Rdel_transcripts.fasta ]; then
curl -o Rdel_transcripts.fasta	https://fungidb.org/common/downloads/Current_Release/RdelemarRA99-880/fasta/data/FungiDB-45_RdelemarRA99-880_AnnotatedTranscripts.fasta
fi

rm -f sra.txt # make sure this is empty

# download data from this study
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP148808
module load aspera

if [ -d /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208621 ]; then
	ln -s /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208621/SRR7208621_1.fastq.gz Rd_24hr_1.fastq.gz
	ln -s /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208621/SRR7208621_2.fastq.gz Rd_24hr_2.fastq.gz
elif [ ! -f SRR7208621/SRR7208621_1.fastq.gz ]; then
	# this dataset is one of several replicates which are for the 24hr time point of spore germination
	echo "SRR7208621" > sra.txt
	/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
	ln -s SRR7208621/SRR7208621_1.fastq.gz Rd_24hr_1.fastq.gz
	ln -s SRR7208621/SRR7208621_2.fastq.gz Rd_24hr_2.fastq.gz
fi

if [ -d /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208587 ]; then
	ln -s /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208587/SRR7208587_1.fastq.gz Rd_3hr_1.fastq.gz
	ln -s /bigdata/stajichlab/shared/tmp/Rhizopus/SRR7208587/SRR7208587_2.fastq.gz Rd_3hr_2.fastq.gz
elif [ ! -f SRR7208587/SRR7208587_1.fastq.gz ]; then
# this dataset is one of several replicates which are for the 3hr time point of spore germination
	echo "SRR7208587" > sra.txt
	/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
	# link the SRA accession to a more useful name
	ln -s SRR7208587/SRR7208587_1.fastq.gz Rd_3hr_1.fastq.gz
	ln -s SRR7208587/SRR7208587_1.fastq.gz Rd_3hr_2.fastq.gz
fi
