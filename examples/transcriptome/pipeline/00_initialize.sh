#!/usr/bin/bash

if [ ! -f Rdel_transcripts.fasta ]; then
	curl -o Rdel_transcripts.fasta	https://fungidb.org/common/downloads/Current_Release/RdelemarRA99-880/fasta/data/FungiDB-45_RdelemarRA99-880_AnnotatedTranscripts.fasta
fi

# download data from this study
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP148808
FASTQ=input
mkdir -p $FASTQ
ln -s /bigdata/stajichlab/shared/tmp/Rhizopus/SRR*/*.fastq.gz $FASTQ/
