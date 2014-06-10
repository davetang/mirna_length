#!/bin/bash

#download mature miRNAs
#Mon Jun  9 14:33:00 JST 2014
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

zcat mature.fa.gz | perl -nle 'if (/^>/){print} else { s/U/T/g; print }' | gzip > mature_thymine.fa.gz

for i in {15..30}
   do seq_by_markov_chain.R $i
   seq_by_equal.R $i
	seq_by_gen_freq.R $i
   random_seq.pl 1000000 $i
done

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip hg38.fa.gz

wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2
tar -xjf bwa-0.7.9a.tar.bz2
cd bwa-0.7.9a/
make

bwa index hg38.fa

for file in `ls my_random*.fa`;
   do echo $file;
   base=`basename $file .fa`
   bwa aln -t 12 hg38.fa $file > $base.sai
   bwa samse hg38.fa $base.sai $file > $base.sam
   samtools view -bS $base.sam > $base.bam
   samtools sort $base.bam ${base}_sorted
   rm $base.sam $base.bam $base.sai
done

#summaries
#perfectly mapped sequences for equal probability of nucleotides
for file in `ls my_random_seq_equal_*sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#perfectly mapped sequences for genome matched probabilities of nucleotides
for file in `ls my_random_seq_gen_freq_*_sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#perfectly mapped sequences for Markov chain generated sequences
for file in `ls my_random_seq_[1-3][0-9]_sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#mapped sequences for equal probability of nucleotides
for file in `ls my_random_seq_equal_*sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done

#mapped sequences for genome matched probabilities of nucleotides
for file in `ls my_random_seq_gen_freq_*_sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done

#mapped sequences for Markov chain generated sequences
for file in `ls my_random_seq_[1-3][0-9]_sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done

