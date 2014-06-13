#!/bin/bash

usage(){
   echo "Usage: $0 <number to generate>"
   exit 1
}

[[ $# -eq 0 ]] && usage
my_number=$1

#download mature miRNAs
if [ ! -f mature.fa.gz ]
then
   wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
   zcat mature.fa.gz | perl -nle 'if (/^>/){print} else { s/U/T/g; print }' > mature_thymine.fa
fi

#download BWA
if [ ! -f bwa*.tar.bz2 ]
then
   wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.9a.tar.bz2
   tar -xjf bwa-0.7.9a.tar.bz2
   cd bwa-0.7.9a/
   make
   cd ..
fi

#download twoBitToFa
if [ ! -f twoBitToFa ]
then
   wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
   chmod 755 twoBitToFa
fi

human=hg38
mouse=mm10
celegans=ce10
zebrafish=danRer7

#create genome directories
if [ ! -d $human ]
then
   mkdir $human
fi
if [ ! -d $mouse ]
then
   mkdir $mouse
fi
if [ ! -d $zebrafish ]
then
   mkdir $zebrafish
fi
if [ ! -d $celegans ]
then
   mkdir $celegans
fi

#download genomes
if [ ! -f $human/$human.2bit ]
then
   wget http://hgdownload.cse.ucsc.edu/goldenPath/$human/bigZips/$human.2bit -O $human/$human.2bit
   twoBitToFa $human/$human.2bit $human/$human.fa
fi
if [ ! -f $mouse/$mouse.2bit ]
then
   wget http://hgdownload.cse.ucsc.edu/goldenPath/$mouse/bigZips/$mouse.2bit -O $mouse/$mouse.2bit
   twoBitToFa $mouse/$mouse.2bit $mouse/$mouse.fa
fi
if [ ! -f $zebrafish/$zebrafish.2bit ]
then
   wget http://hgdownload.cse.ucsc.edu/goldenPath/$zebrafish/bigZips/$zebrafish.2bit -O $zebrafish/$zebrafish.2bit
   twoBitToFa $zebrafish/$zebrafish.2bit $zebrafish/$zebrafish.fa
fi
if [ ! -f $celegans/$celegans.2bit ]
then
   wget http://hgdownload.cse.ucsc.edu/goldenPath/$celegans/bigZips/$celegans.2bit -O $celegans/$celegans.2bit
   twoBitToFa $celegans/$celegans.2bit $celegans/$celegans.fa
fi

#index genomes if all bwa index files are missing
if [ ! -f $zebrafish/$zebrafish.fa.amb ] && [ ! -f $zebrafish/$zebrafish.fa.ann ] && [ ! -f $zebrafish/$zebrafish.fa.bwt ] && [ ! -f $zebrafish/$zebrafish.fa.pac ] && [ ! -f $zebrafish/$zebrafish.fa.sa ]
then
   bwa-0.7.9a/bwa index $zebrafish/$zebrafish.fa
fi

if [ ! -f $mouse/$mouse.fa.amb ] && [ ! -f $mouse/$mouse.fa.ann ] && [ ! -f $mouse/$mouse.fa.bwt ] && [ ! -f $mouse/$mouse.fa.pac ] && [ ! -f $mouse/$mouse.fa.sa ]
then
   bwa-0.7.9a/bwa index $mouse/$mouse.fa
fi

if [ ! -f $human/$human.fa.amb ] && [ ! -f $human/$human.fa.ann ] && [ ! -f $human/$human.fa.bwt ] && [ ! -f $human/$human.fa.pac ] && [ ! -f $human/$human.fa.sa ]
then
   bwa-0.7.9a/bwa index $human/$human.fa
fi

if [ ! -f $celegans/$celegans.fa.amb ] && [ ! -f $celegans/$celegans.fa.ann ] && [ ! -f $celegans/$celegans.fa.bwt ] && [ ! -f $celegans/$celegans.fa.pac ] && [ ! -f $celegans/$celegans.fa.sa ]
then
   bwa-0.7.9a/bwa index $celegans/$celegans.fa
fi

if [ ! -f hg38_init_prob.Robject ] || [ ! -f hg38_trans_mat.Robject ] || [ ! -f mm10_init_prob.Robject ] || [ ! -f mm10_trans_mat.Robject ] || [ ! -f ce10_init_prob.Robject ] || [ ! -f ce10_trans_mat.Robject ] || [ ! -f danRer7_init_prob.Robject ] || [ ! -f danRer7_trans_mat.Robject ]
then
   R --no-save < dinucleotide.R
fi

if [ ! -f hg38_nuc_freq.Robject ] || [ ! -f ce10_nuc_freq.Robject ] || [ ! -f danRer7_nuc_freq.Robject ] || [ ! -f mm10_nuc_freq.Robject ]
then
   R --no-save < genome_freq.R
fi

#generate random sequences
for org in $human $mouse $celegans $zebrafish
   do for i in {15..30}
      do
      if [ ! -f $org/my_random_seq_markov_${i}_${my_number}.fa ]
      then
         touch $org/my_random_seq_markov_${i}_${my_number}.fa
         Rscript seq_by_markov_chain.R $i ${org}_trans_mat.Robject ${org}_init_prob.Robject $my_number $org
      fi
      if [ ! -f $org/my_random_seq_equal_${i}_${my_number}.fa ]
      then
         touch $org/my_random_seq_equal_${i}_${my_number}.fa
         Rscript seq_by_equal.R $i $my_number $org
      fi
      if [ ! -f $org/my_random_seq_gen_freq_${i}_${my_number}.fa ]
      then
         touch $org/my_random_seq_gen_freq_${i}_${my_number}.fa
         Rscript seq_by_gen_freq.R $i ${org}_nuc_freq.Robject $my_number $org
      fi
   done
done

#align random sequences
for org in $human $mouse $celegans $zebrafish
   do for file in `ls $org/my_random*.fa`;
      do echo $file;
      base=`basename $file .fa`
      if [ ! -f $org/${base}_sorted.bam ]
      then
         touch $org/${base}_sorted.bam
         bwa-0.7.9a/bwa aln -t 12 $org/$org.fa $file > $org/$base.sai
         bwa-0.7.9a/bwa samse $org/$org.fa $org/$base.sai $file > $org/$base.sam
         samtools view -bS $org/$base.sam > $org/$base.bam
         samtools sort $org/$base.bam $org/${base}_sorted
         rm $org/$base.sam $org/$base.bam $org/$base.sai
      fi
   done
done

#transition images
if [ ! -f image/$celegans.eps ]
then
   Rscript image/transition.R ${ce10}_trans_mat.Robject
   convert image/$celegans.eps -trim +repage image/$celegans.pdf
fi
if [ ! -f image/$zebrafish.eps ]
then
   Rscript image/transition.R ${zebrafish}_trans_mat.Robject
   convert image/$zebrafish.eps -trim +repage image/$zebrafish.pdf
fi
if [ ! -f image/$human.eps ]
then
   Rscript image/transition.R ${human}_trans_mat.Robject
   convert image/$human.eps -trim +repage image/$human.pdf
fi
if [ ! -f image/$human.eps ]
then
   Rscript image/transition.R ${mouse}_trans_mat.Robject
   convert image/$mouse.eps -trim +repage image/$mouse.pdf
fi

#miRNA lengths
if [ -f image/celegans_mirbase_length.eps ]
then
   convert image/celegans_mirbase_length.eps image/celegans_mirbase_length.pdf
fi
if [ -f image/mouse_mirbase_length.eps ]
then
   convert image/mouse_mirbase_length.eps image/mouse_mirbase_length.pdf
fi
if [ -f image/human_mirbase_length.eps ]
then
   convert image/human_mirbase_length.eps image/human_mirbase_length.pdf
fi
if [ -f image/zebrafish_mirbase_length.eps ]
then
   convert image/zebrafish_mirbase_length.eps image/zebrafish_mirbase_length.pdf
fi

#perfectly mapped
if [ -f image/ce10_perfect_mapped_1000000.eps ]
then
   convert image/ce10_perfect_mapped_1000000.eps image/ce10_perfect_mapped_1000000.pdf
fi
if [ -f image/danRer7_perfect_mapped_1000000.eps ]
then
   convert image/danRer7_perfect_mapped_1000000.eps image/danRer7_perfect_mapped_1000000.pdf
fi
if [ -f image/hg38_perfect_mapped_1000000.eps ]
then
   convert image/hg38_perfect_mapped_1000000.eps image/hg38_perfect_mapped_1000000.pdf
fi
if [ -f image/mm10_perfect_mapped_1000000.eps ]
then
   convert image/mm10_perfect_mapped_1000000.eps image/mm10_perfect_mapped_1000000.pdf
fi

#mapped
if [ -f image/ce10_mapped_1000000.eps ]
then
   convert image/ce10_mapped_1000000.eps image/ce10_mapped_1000000.pdf
fi
if [ -f image/danRer7_mapped_1000000.eps ]
then
   convert image/danRer7_mapped_1000000.eps image/danRer7_mapped_1000000.pdf
fi
if [ -f image/hg38_mapped_1000000.eps ]
then
   convert image/hg38_mapped_1000000.eps image/hg38_mapped_1000000.pdf
fi
if [ -f image/mm10_mapped_1000000.eps ]
then
   convert image/mm10_mapped_1000000.eps image/mm10_mapped_1000000.pdf
fi

exit
