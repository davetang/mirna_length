#!/bin/bash

human=hg38

#align random sequences
for file in `ls $human/my_random*.bam`;
   do echo $file
   exon=`intersectBed -abam $file -b ../defining_genomic_regions/gencode_v19_exon_merged.bed.gz -bed -wo | cut -f16 | perl -nle '$n+=$_; END { print $n }'`
   intron=`intersectBed -abam $file -b ../defining_genomic_regions/gencode_v19_intron.bed.gz -bed -wo | cut -f16 | perl -nle '$n+=$_; END { print $n }'`
   intergenic=`intersectBed -abam $file -b ../defining_genomic_regions/gencode_v19_intergenic.bed.gz -bed -wo | cut -f16 | perl -nle '$n+=$_; END { print $n }'`
   total=`expr $exon + $intron + $intergenic`
   exon_percent=`bc -l<<<$exon*100/$total`
   intron_percent=`bc -l<<<$intron*100/$total`
   intergenic_percent=`bc -l<<<$intergenic*100/$total`
   printf "Exon %.2f\n" $exon_percent
   printf "Intron %.2f\n" $intron_percent
   printf "Intergenic %.2f\n" $intergenic_percent
done

echo Done

exit
