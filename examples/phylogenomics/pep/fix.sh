
for file in $(ls *.proteins.fa)
do
	base=$(basename $file .proteins.fa)
	perl -p -e 's/^>([^_]+)_/>$1|$1_/' $file > $base.aa.fasta
done
