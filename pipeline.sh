#!/bin/bash

#download mature miRNAs
#Mon Jun  9 14:33:00 JST 2014
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

zcat mature.fa.gz | perl -nle 'if (/^>/){print} else { s/U/T/g; print }' | gzip > mature_thymine.fa.gz

#download genomes
#human
mkdir hg38
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O human/hg38.fa.gz
gunzip human/hg38.fa.gz
#mouse
mkdir mm10
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit -O mouse/mm10.2bit
twoBitToFa mouse/mm10.2bit mouse/mm10.fa
#zebrafish
mkdir danRer7
wget http://hgdownload.cse.ucsc.edu/goldenPath/danRer7/bigZips/danRer7.fa.gz -O zebrafish/danRer7.fa.gz
gunzip zebrafish/danRer7.fa.gz
#celegans
mkdir ce10
wget http://hgdownload.cse.ucsc.edu/goldenPath/ce10/bigZips/ce10.2bit -O celegans/ce10.2bit
twoBitToFa celegans/ce10.2bit celegans/ce10.fa

#download BWA
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2
tar -xjf bwa-0.7.9a.tar.bz2
cd bwa-0.7.9a/
make

#index genomes
bwa index danRer7/danRer7.fa
bwa index ce10/ce10.fa
bwa index hg38/hg38.fa
bwa index mm10/mm10.fa

R --no-save < dinucleotide.R
R --no-save < genome_freq.R

#generate random sequences
for org in hg38 mm10 ce10 danRer7
   do for i in {15..30}
      do seq_by_markov_chain.R $i ${org}_trans_mat.Robject ${org}_init_prob.Robject 20 $org
      seq_by_equal.R $i 20 $org
   	seq_by_gen_freq.R $i ${org}_nuc_freq.Robject 20 $org
   done
done

#align random sequences
for org in hg38 mm10 ce10 danRer7
   do for file in `ls $org/my_random*.fa`;
      do echo $org/$file;
      base=`basename $org/$file .fa`
      bwa aln -t 12 $org/$org.fa $org/$file > $org/$base.sai
      bwa samse $org/$org.fa $org/$base.sai $org/$file > $org/$base.sam
      samtools view -bS $org/$base.sam > $org/$base.bam
      samtools sort $org/$base.bam $org/${base}_sorted
      rm $org/$base.sam $org/$base.bam $org/$base.sai
   done
done

exit

#summaries
#perfectly mapped sequences for equal probability of nucleotides
#for file in `ls my_random_seq_equal_*sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#perfectly mapped sequences for genome matched probabilities of nucleotides
#for file in `ls my_random_seq_gen_freq_*_sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#perfectly mapped sequences for Markov chain generated sequences
#for file in `ls my_random_seq_[1-3][0-9]_sorted.bam`; do echo $file; bam_stat.pl $file | grep "^0"; done

#mapped sequences for equal probability of nucleotides
#for file in `ls my_random_seq_equal_*sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done

#mapped sequences for genome matched probabilities of nucleotides
#for file in `ls my_random_seq_gen_freq_*_sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done

#mapped sequences for Markov chain generated sequences
#for file in `ls my_random_seq_[1-3][0-9]_sorted.bam`; do echo $file; bam_stat.pl $file | grep ^Mapped: | cut -f2 -d' '; done
