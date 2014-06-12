#!/bin/bash

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

#index genomes
if [ ! -f $zebrafish/$zebrafish.fa.amb ] || [ ! -f $zebrafish/$zebrafish.fa.ann ] || [ ! -f $zebrafish/$zebrafish.fa.bwt ] || [ ! -f $zebrafish/$zebrafish.fa.pac ] || [ ! -f $zebrafish/$zebrafish.fa.sa ]
then
   bwa-0.7.9a/bwa index $zebrafish/$zebrafish.fa
fi

if [ ! -f $mouse/$mouse.fa.amb ] || [ ! -f $mouse/$mouse.fa.ann ] || [ ! -f $mouse/$mouse.fa.bwt ] || [ ! -f $mouse/$mouse.fa.pac ] || [ ! -f $mouse/$mouse.fa.sa ]
then
   bwa-0.7.9a/bwa index $mouse/$mouse.fa
fi

if [ ! -f $human/$human.fa.amb ] || [ ! -f $human/$human.fa.ann ] || [ ! -f $human/$human.fa.bwt ] || [ ! -f $human/$human.fa.pac ] || [ ! -f $human/$human.fa.sa ]
then
   bwa-0.7.9a/bwa index $human/$human.fa
fi

if [ ! -f $celegans/$celegans.fa.amb ] || [ ! -f $celegans/$celegans.fa.ann ] || [ ! -f $celegans/$celegans.fa.bwt ] || [ ! -f $celegans/$celegans.fa.pac ] || [ ! -f $celegans/$celegans.fa.sa ]
then
   bwa-0.7.9a/bwa index $celegans/$celegans.fa
fi

R --no-save < dinucleotide.R
R --no-save < genome_freq.R

my_number=1000000

#generate random sequences
for org in $human $mouse $celegans $zebrafish
   do for i in {15..30}
      do Rscript seq_by_markov_chain.R $i ${org}_trans_mat.Robject ${org}_init_prob.Robject $my_number $org
      Rscript seq_by_equal.R $i $my_number $org
      Rscript seq_by_gen_freq.R $i ${org}_nuc_freq.Robject $my_number $org
   done
done

#align random sequences
for org in $human $mouse $celegans $zebrafish
   do for file in `ls $org/my_random*.fa`;
      do echo $org/$file;
      base=`basename $org/$file .fa`
      bwa-0.7.9a/bwa aln -t 12 $org/$org.fa $org/$file > $org/$base.sai
      bwa-0.7.9a/bwa samse $org/$org.fa $org/$base.sai $org/$file > $org/$base.sam
      samtools view -bS $org/$base.sam > $org/$base.bam
      samtools sort $org/$base.bam $org/${base}_sorted
      rm $org/$base.sam $org/$base.bam $org/$base.sai
   done
done

#transition images
Rscript image/transition.R ce10_trans_mat.Robject
Rscript image/transition.R danRer7_trans_mat.Robject
Rscript image/transition.R hg38_trans_mat.Robject
Rscript image/transition.R mm10_trans_mat.Robject
convert image/ce10.eps -trim +repage image/ce10.pdf
convert image/danRer7.eps -trim +repage image/danRer7.pdf
convert image/hg38.eps -trim +repage image/hg38.pdf
convert image/mm10.eps -trim +repage image/mm10.pdf

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
