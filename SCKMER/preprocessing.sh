#!/bin/bash

FASTADIR=barcode_fasta
COUNTDIR=barcode_count_30mer

FINALBARCODELABEL=30mer.barcode.label.selected.txt
FINALCOUNTTABLE=30mer.count.table.selected.npz

mkdir ${COUNTDIR}
mv ${FASTADIR}/*.kmer.txt ${COUNTDIR}/

for i in `ls ${FASTADIR}/*.fasta`;do n=`basename $i`; n=${n/.fasta/}; l=`grep -c \> $i`; echo $n $l;done > n.seq.txt

cd ${COUNTDIR}

wc -l *.txt  | grep -v total| perl -pe 's/^\s+(\d+)\s+(\S+?)\.kmer.txt$/$2\t$1/' > ../n.kmer.txt

cd ..

echo -e "barcode\tn.kmers\tn.seqs" > n.kmers.seqs.txt

perl -e 'open(K,"<$ARGV[0]");open(S,"<$ARGV[1]");while(<K>){chomp;@_=split;$k{$_[0]}=$_[1]}while(<S>){chomp;@_=split;$s{$_[0]}=$_[1]}close(K);close(S);foreach (keys%k){print $_,"\t",$k{$_},"\t",$s{$_},"\n"}' n.kmer.txt n.seq.txt >> n.kmers.seqs.txt

Rscript plot.R

perl -ne 'chomp;@_=split;if($_[2]>=20000 && $_[2]<=150000){ print $_,"\n"}' n.kmers.seqs.txt > n.kmers.seqs.kept.txt

rm -f uniq.kmer.tmp

touch uniq.kmer.tmp

cut -f 1 n.kmers.seqs.kept.txt| while read bc; do
    cut -f 1 ${COUNTDIR}/${bc}.kmer.txt >> uniq.kmer.tmp;
done

echo "kmer" > uniq.kmer.txt

sort -u uniq.kmer.tmp >> uniq.kmer.txt

rm -f uniq.kmer.tmp

./remove_simplerepeat.py > uniq.kmer.filtered.txt

./build_sp.py ${COUNTDIR}

mv barcode.label.txt ${FINALBARCODELABEL}
mv count.table.npz ${FINALCOUNTTABLE}
