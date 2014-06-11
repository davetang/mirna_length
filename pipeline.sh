#!/bin/bash

#download mature miRNAs
#Mon Jun  9 14:33:00 JST 2014
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

zcat mature.fa.gz | perl -nle 'if (/^>/){print} else { s/U/T/g; print }' | gzip > mature_thymine.fa.gz

#download genomes
#human
mkdir human
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O human/hg38.fa.gz
gunzip human/hg38.fa.gz
#mouse
mkdir mouse
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit -O mouse/mm10.2bit
twoBitToFa mouse/mm10.2bit mouse/mm10.fa
#zebrafish
mkdir zebrafish
wget http://hgdownload.cse.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.gz -O zebrafish/danRer7.fa.gz
gunzip zebrafish/danRer7.fa.gz
#celegans
mkdir celegans
wget http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit -O celegans/ce10.2bit
twoBitToFa celegans/ce10.2bit celegans/ce10.fa

#download BWA
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2
tar -xjf bwa-0.7.9a.tar.bz2
cd bwa-0.7.9a/
make

#index genomes
bwa index zebrafish/danRer7.fa
bwa index celegans/ce10.fa
bwa index human/hg38.fa
bwa index mouse/mm10.fa

#generate random sequences
for i in {15..30}
   do seq_by_markov_chain.R $i
   seq_by_equal.R $i
	seq_by_gen_freq.R $i
done

#align random sequences
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

