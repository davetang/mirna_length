#!/bin/bash

usage(){
   echo "Usage: $0 <number to generate>"
   exit 1
}

[[ $# -eq 0 ]] && usage
my_number=$1
my_length=22

human=hg38
org=$human

#generate random sequences
for i in a b c d e
   do
   if [ ! -f my_random_seq_equal_${my_length}_${my_number}_$i.fa ]
   then
      touch my_random_seq_equal_${my_length}_${my_number}_$i.fa
      Rscript ../script/seq_by_equal.R $my_length $my_number .
      mv my_random_seq_equal_${my_length}_${my_number}.fa my_random_seq_equal_${my_length}_${my_number}_$i.fa
   fi

   if [ ! -f my_random_seq_gen_freq_${my_length}_${my_number}_$i.fa ]
   then
      touch my_random_seq_gen_freq_${my_length}_${my_number}_$i.fa
      Rscript ../script/seq_by_gen_freq.R $my_length ../data/${org}_nuc_freq.Robject $my_number .
      mv my_random_seq_gen_freq_${my_length}_${my_number}.fa my_random_seq_gen_freq_${my_length}_${my_number}_$i.fa
   fi

   if [ ! -f my_random_seq_markov_${my_length}_${my_number}_$i.fa ]
   then
      touch my_random_seq_markov_${my_length}_${my_number}_$i.fa
      Rscript ../script/seq_by_markov_chain.R $my_length ../data/${org}_trans_mat.Robject ../data/${org}_init_prob.Robject $my_number .
      mv my_random_seq_markov_${my_length}_${my_number}.fa my_random_seq_markov_${my_length}_${my_number}_$i.fa
   fi
done

exit

#align random sequences
do for file in `ls my_random*.fa`;
   do base=`basename $file .fa`
   ../bwa-0.7.9a/bwa aln -t 12 ../$org/$org.fa $file > $base.sai
   ../bwa-0.7.9a/bwa samse ../$org/$org.fa $base.sai $file > $base.sam
   samtools view -bS $base.sam > $base.bam
   samtools sort $base.bam ${base}_sorted
   rm $base.sam $base.bam $base.sai
done

echo Done

exit
