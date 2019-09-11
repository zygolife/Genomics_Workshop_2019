# Retrieving Sequences from databases

Some of this is mentioned in the [Genome Assembly](Assembly) section. But more detail about different ways to retrieve sequences is locate here.

# Remote databases

The International Nucleotide Sequence Database Collaboration ([INSDC](http://www.insdc.org/)) are the joint databases for sequence and related biomedical data deposition. These are critical central tools for archive of sequence data for the scientific community.  Often when we say "deposited in GenBank" we mean this central repository which, but there are data for Sequence Reads, Assemblies, annotated genomes, Gene Expression, and Individual sequence records.

## Downloading FASTA databases from NCBI, Uniprot

On HPCC there are already databases installed and indexed for BLAST searches.

```BASH
#SBATCH -N 1 -n 16
module load db-ncbi
module load ncbi-blast
# loads the current ncbi folder as env variables
# $BLASTDB and $NCBI_DB
# after loading this you can run blast without specifying
# the location of the databases
blastp -db nr -query seqs.fasta -out seqs-nr.blastp -num_threads 16 -evalue 1e-5
```

### FTP / Web downloads

**NCBI**

The NCBI databases include several useful resources. Note these are now giant databases in some places so care is needed in whether you can download these to your own folder and if it already exists on cluster you should try to use those.

* [swissprot](curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/swissprot.gz)
* [nr](curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/nr.gz) - non-redundant protein seqs, very large ...
* [nt](curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/nt.gz) - non-redundant nucleotide seqs, very large ...
* [env_nr](curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/env_nr.gz) - enviromental protein seqs
* [env_nt](curl -O ftp://ftp.ncbi.nih.gov:/blast/db/FASTA/env_nt.gz) - enviromental nucl seqs

Another resource that is helpful is [Refseq](ftp://ftp.ncbi.nih.gov/refseq/release) which are _somewhat verified_ sequences from genomes.

To for example get all fungal refseq proteins use the lftp tool. Here is an interactive session:
```BASH
lftp ftp://ftp.ncbi.nih.gov/refseq/release
lftp> cd fungi
lftp> mget fungi.*.protein.faa.gz
lftp> exit
pigz -dc *.faa.gz > refseq_fungi.faa
```

**Uniprot**

To Download uniprot_swissprot database via ftp protocol. See https://www.uniprot.org/downloads for more download files available including uniref which is a set of clustered sequences at 100%, 90%, and 50% identity which can reduce the size of the total protein database but still leave representatives.

* [uniprot_swissprot](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
* [uniref50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz)

The [UniRef50 database](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/) for this as it isn't too big but useful for some relatively fast searching and more comprehensive than swissprot for taxonomic representation.

## Downloading from SRA

The easiest way to download from SRA is using the [sra_download.pl](https://github.com/gerthmicha/perlscripts/blob/master/sra_download.pl) script. This is already installed on the cluster in `/bigdata/stajichlab/shared/bin/sra_download.pl`.

To run this create a file that has a line for each SRA accession number. I have called it `sra.txt` here. This will download the fastq for all the data sets in the file and create a folder for each one.
```bash
#SBATCH -p short -N 1 -n 2 --mem 4gb
module load aspera
/bigdata/stajichlab/shared/bin/sra_download.pl --ascp --id $ASPERAKEY sra.txt
```

## Downloading sequence records from GenBank

You can use several tools to download accessions from genbank.  It does require certain versions of perl are installed or conda.


```BASH
module load perl/5.20.2
bp_download_query_genbank.pl --query 'AY295118.1'
>AY295118 Parmelia ernstiae voucher MAF 9805 tubulin gene, partial cds.
GAGGACATTCCTCCATAATGTGATACGTAGCTCACAGCTTTCAAGGCTTCAAACAACAAA
TATGTTCCTCGTGCCGTACTCGTCGATCTCGAGCCTGGTACCATGGATGCTGTCCGCGCT
GGTCCTTTTGGCCAGCTTTTCCGACCCGATAACTTCGTATTTGGTCAATCTGGTGCTGGT
AATAATTGGGCTAAGGGTCATTACACCGAGGGTGCAGAATTGGTGGACCAAGTCCTCGAT
GTTGTGCGTCGAGAGGCTGAAGGATGCGACTGCCTCCAGGGCTTCCAGATCACGCACTCC
CTCGGTGGTGGAACTGGTGCTGGTATGGGTACGCTTTTGATCTCGAAAATCCGTGAGGAG
TTCCCAGATCGTATGATGGCTACATTCTCCGTGGTTCCTTCACCAAAGGTATCCGACACT
GTTGTGGAGCCATACAACGCTACTCTCTCCGTGCATCAATTGGTCGAGAACTCGGATGAG
ACCTTCTGTATCGATAATGAGGTTGGTCAAGTGCGATTTTTTCACAGAGGCGCAAGGACT
GATATGTCAATCTAGGCGCTCTATGACATTTGCATGCGCACCCTCAAGCTCTCCAACCCA
TCCTACGGGGATCTTAACCACCTTGTCTCCGCGGTCATGTCTGGTGTTACCACCTGCCTC
CGTTTCCCCGGTCAACTCAATTCCGACCTTCGAAAACTAGCCGTCAACATGGTCCCATTT
CCCCGTCTACATTTCTTCATGGTTGGCTTCGCACCTCTTACCAGCCGAGGTGCTAACTCA
TTCCGTGCGGTCAGCGTACCAGAATTGACCCAACAAATGTACGAC
```

If you want to retrieve a number of sequences at a time you can specify a query. Below are the options for running the tool. If you want to retrieve data from protein database you need to specify the database with `--db` option.

```
bp_download_query_genbank --query "Neurospora[ORGN]" --db nucest -o Ncrassa_ESTs.fa --format fasta


Other options
 Provide ONE of:

  -q --query query string OR
  --queryfile profile file with query OR
  --gi --gis --gifile file with list of GIs to download

 Database type:

 -d --db database (nucleotide [default], nucest, protein, )

 -o --out --outfile output file (results are displayed on screen otherwise)
 -f --format sequence file output format (fasta by default)
 -v --verbose debugging output

Query options
 --maxids maximum number of IDs to retrieve in a set (100 at a time by default)
 --reldate
 --maxdate maxdate for a record
 --mindate minimum date for record
 --datetype edat or mdat (entered or modified)
```

## Specialized fungal databases

### FungiDB

The [FungiDB](https://fungidb.org/fungidb/) project provides access to a set of Fungal genomes loaded into this system. The resources for downloads are available at [this link](https://fungidb.org/common/downloads/) which includes current and previous releases.  Data sets are organized by Abbreviations of genus + species and strain name.  For example the Genome, CDS, Protein, and Transcripts associated with the _Neurospora crassa_ OR74A strain are available from [this link](https://fungidb.org/common/downloads/Current_Release/NcrassaOR74A/fasta/data/).

### JGI

The JGI [Mycocosm](http://mycocosm.jgi.doe.gov) provides one of the largest collection of fungal genomes through the sequencing and annotation project. There are more than [1000 genomes available](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=fungi) through [several interfaces](https://mycocosm.jgi.doe.gov/fungi/fungi.info.html) hosted by the [JGI](http://jgi.doe.gov). [Some scripts](https://github.com/1KFG/2019_dataset) for automation of downloads are needed to directly extract data from the site onto linux clusters. GLOBUS and other functionality do exist for dataset downloads as well.

### Ensembl

Ensembl provides a nearly complete set of public deposited genomes into GenBank organized by major domains. The [Ensembl Fungi](http://fungi.ensembl.org/index.html) is a portal with access to thousands

### Saccharomyces or Candida Genome Database

See [SGD](http://yeastgenome.org) and [CGD](http://candidagenome.org) for main site for genome browsers, comparative tools, and access to primary sequence data associated with these fungi.

The FTP site for yeast see ftp://ftp.yeastgenome.org/sequence/S288C_reference for example which has access to the yeast ORFs [proteins](ftp://ftp.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz) and [coding sequence](ftp://ftp.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz) as well as many other resources like upstream promotor files and other features. Data from multiple _Candida_ species is available from http://candidagenome.org/download/sequence/

# Local databases

Once you have data files downloaded you can use these.

## FASTA Files

### HMMER esl-sfetch

I find esl-sfetch one of the better tools for fasta file indexing. Database must be indexed first.

```BASH
module load hmmer
esl-sfetch --index database.fasta
# fetch record based on single ID passed in on cmdline
esl-sfetch database.fasta accession > accession.fa
# fetch multiple records based on a list of IDs passed in a file (note the -f option)
esl-sfetch -f database.fasta list_of_ids > seqs.fa
# fetch list of IDs passed in on STDIN using a pipe and specifying the input file as '-'
cat ids_to_fetch | esl-sfetch -f database.fasta - > seqs.fa
```

### cdbfasta
cdbfasta ([Constant database](https://github.com/gpertea/cdbfasta)) is a useful for indexing fasta and fastq files for retrieval by sequence ID.

```BASH
module load cdbfasta
cdbfasta database.fasta

# to retrieve
echo "A" | cdbyank databasa.fasta.cidx  > A.fa
cat list_of_ids | cdbyank database.fasta.cidx > retrieved.fa
```

### samtools (1.9 and later)

Samtools provides indexing and retrieval of FASTA

```bash
#SBATCH -p short -N 1 -n 4
module load samtools

# index file
samtools faidx DNA_sequences.fasta
```

To retrieve a sequence read after the file is indexed, where accession is the first text after the > in FASTA file, eg `scaffold_1` is the accession in the following:
```text
>scaffold_1
TGCATGTCTAAGTATAAGCAATTATACCGTGAAACTGCGAATGGCTCATTAAATCAGTTATCGTTTATTTGATAGTACCTTACTACTTGGATAACCGTGGTAATTCTAGAGCTAATACATGCTGAAAACCCCAACTTCGGGAGGGGTGTATTTATTAGATAAAAAACCAACGCCCTTCGGGGCTTCTTGGTGATTCATGATAACTTTACGGATCGCATGGCCTTGCGCCGGCGACGGTTCATTCAAATTTCTGCCCTATCAACTTTCGATGGTAAGGTATTGGCTTACCATGGTTTCAACGGGTAACGGGGAATTAGGGTTCGATTCCGGAGAGGGAGCCTGAGAAACGGCTACCACA
```
To retrieve this sequence from the indexed file use this (specify DNA_sequences.fasta if you used the uncompressed file).
```
module load samtools
samtools faidx DNA_sequences.fasta scaffold_1
```
### BLAST indexing

The BLAST indexing to setup a database for sequence alignment and searching also allows retrieval of sequences by identifier.

```BASH
module load ncbi-blast
# index a nucleotide database
# to index a protein database change -dbtype from 'nucl' to 'prot'
makeblastdb -in sequences.fasta -dbtype nucl -parse_seqids

# to retrieve sequences
blastdbcmd -entry ACCESSION -db sequences.fasta -out ACCESSION.fasta

# use this database for sequence searches
# report the output as tab delimited format (outfmt 6)
blastn -query myquery.fasta -db sequences.fasta -out myquery-vs-seqs.BLASTN -outfmt 6 -evalue 1e-5

# Do a protein db search
blastp -query myquery.fasta -db protseqdb.fasta -out myquery-vs-seqs.BLASTP -outfmt 6 -evalue 1e-4

# many other options for BLAST using blastx, tblastn, tblastx and many more options for running BLAST not explained here.
```


### DIAMOND indexing

[DIAMOND](https://github.com/bbuchfink/diamond) is a rapid aligner for protein and translated searches which can operate on short sequence reads as well as assembled genomes.  

DIAMOND does not provide a way to extract sequences back out from these indexed databases. Will report a tab delimited output file (m8 stands for the OLD NCBI -mformat 8 output which is tab delimited).

```BASH
module load diamond
makedb --in my_protein_db.fasta -d mydb
diamond blastx -d mydb -q reads.fna -o hits.m8
```

### Short read aligner database indexing

Indexing DNA database for aligning short read DNA sequences against this database (usually a genome).

Indexing for *bwa* in order to setup searches:

```bash
module load bwa
bwa index database.fa
```
Indexing for *bowtie2*:

```bash
module load bowtie2
bowtie2-build database.fa database
```

Indexing for *gmap/gsnap*:

```BASH
module load gmap
gmap_build -D genome_index -d genome_name database.fa
```

Indexing for *kallisto* (RNASeq analysis):

```BASH
module load kallisto
kallisto index -i transcripts.idx transcripts.fasta
```

## FASTQ Files

### cdbfasta

cdbfasta ([Constant database](https://github.com/gpertea/cdbfasta)) is a useful for indexing fasta and fastq files for retrieval by sequence ID.

```BASH
module load cdbfasta
cdbfasta -Q reads.fastq

# retrieve seqs by
echo "ACCESSION" | cdbyank reads.fastq.cidx > fetched_read.fq
cat list_of_ids | cdbyank reads.fastq.cidx > retrieved.fq
```

### samtools

Samtools provides indexing and retrieval of FASTQ Files.

If the file is compressed (.gz) it must be compressed with the bgzip tool - which is part of the htslib package. So if the file exists already as a compressed file you need to uncompress and recompress with bgzip.

```bash
#SBATCH -p short -N 1 -n 4
module load samtools
pigz -d READFILE.fq.gz
bgzip --threads 4 READFILE.fq

# now index
samtools fqidx READFILE.fq.gz

# you can also index an uncompressed file
samtools fqidx READFILE.fq
```
To retrieve a sequence read after the file is indexed, where accession is the first text after the @ in FASTQ file, eg `ERR1309286.4` is the accession in the following:
```text
@ERR1309286.4 H4:C3F32ACXX:2:1101:1849:2436/1
CTCTATTTCATCACGTTCGAGAAGATCGCTACGCTTATCGAATTCCAGATTATCATTGTCCGCTTCAACTTCTAGAGAAACTGTGCATGATAATGAGATGC
+
@CCFFFFFGHHGHJJIIJIHJIJJJJIIIIIJIIFJIIJJGIGIEHGIHIIGJIIIJJJJJJIHGF:BDBEEEEEDEA>>CDDCDDDEDDDEDDDDD<>@C
```
To retrieve this sequence from the indexed file use
```
module load samtools
samtools fqidx READFILE.fq.gz ERR1309286.4
```
