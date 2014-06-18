#!/bin/bash

usage(){
   echo "Usage: $0 <number of sequences generated>"
   exit 1
}

[[ $# -eq 0 ]] && usage
my_number=$1

human=hg38
mouse=mm10
celegans=ce10
zebrafish=danRer7

#mapping summaries for perfectly mapped
#equal probability of nucleotides
for org in $human $mouse $celegans $zebrafish
   do echo $org
   rm -f $org/*$my_number.txt
   for i in {15..30}
      do echo $i
      #count equal
      if [ -f $org/my_random_seq_equal_${i}_${my_number}_sorted.bam ]
      then
         script/bam_stat.pl $org/my_random_seq_equal_${i}_${my_number}_sorted.bam | grep "^0" | cut -f2 >> $org/equal_perfect_mapped_$my_number.txt
         script/bam_stat.pl $org/my_random_seq_equal_${i}_${my_number}_sorted.bam | grep ^Mapped: | cut -f2 -d' ' >> $org/equal_mapped_$my_number.txt
      fi
      #count genome freq
      if [ -f $org/my_random_seq_gen_freq_${i}_${my_number}_sorted.bam ]
      then
         script/bam_stat.pl $org/my_random_seq_gen_freq_${i}_${my_number}_sorted.bam | grep "^0" | cut -f2 >> $org/gen_freq_perfect_mapped_$my_number.txt
         script/bam_stat.pl $org/my_random_seq_gen_freq_${i}_${my_number}_sorted.bam | grep ^Mapped: | cut -f2 -d' ' >> $org/gen_freq_mapped_$my_number.txt
      fi
      #count Markov
      if [ -f $org/my_random_seq_markov_${i}_${my_number}_sorted.bam ]
      then
         script/bam_stat.pl $org/my_random_seq_markov_${i}_${my_number}_sorted.bam | grep "^0" | cut -f2 >> $org/markov_perfect_mapped_$my_number.txt
         script/bam_stat.pl $org/my_random_seq_markov_${i}_${my_number}_sorted.bam | grep ^Mapped: | cut -f2 -d' ' >> $org/markov_mapped_$my_number.txt
      fi
      done
done

exit
