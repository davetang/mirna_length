#!/bin/env perl
 
use strict;
use warnings;
 
my $usage = "Usage: $0 <infile.bam>\n";
my $infile = shift or die $usage;
 
my %ed = ();
my %multimappers = ();
my $mapped = '0';
my $unmapped = '0';

if ($infile =~ /\.sam$/){
   open (IN,'<',$infile) || die "Could not open $infile: $!\n";
} else {
   open(IN,'-|',"samtools view $infile") || die "Could not open $infile: $!\n";
}

while(<IN>){
   chomp;
   next if (/^@/);
   my ($junk,$flag) = split(/\t/);

   if ($flag & 0x4){
      ++$unmapped;
      next;
   } else {
      ++$mapped;
	}

   #965850_13       16      chr1    204577  0       13M     *       0       0       GGCTGCCAGAAGA   *       XT:A:N  NM:i:0  XN:i:13 X0:i:582        XM:i:0  XO:i:0  XG:i:0  MD:Z:13
   if (/NM:i:(\d+)/){
      my $ed = $1;
      if (exists $ed{$ed}){
         $ed{$ed}++;
      } else {
         $ed{$ed} = '1';
      }
   } else {
      die "Could not get mismtach stat on line $.: $_\n";
   }
 
   if (/X0:i:(\d+)/){
      my $mapped_time = $1;
      if (exists $multimappers{$mapped_time}){
         $multimappers{$mapped_time}++;
      } else {
         $multimappers{$mapped_time} = '1';
      }
   } else {
      die "Could not get multi mapped stat; X0 tag missing\n";
   }
}
close(IN);
 
print "Mapped: $mapped\n";
print "Unmapped: $unmapped\n";
 
print "Edit distances\n";
foreach my $ed (sort {$a <=> $b} keys %ed){
   print "$ed\t$ed{$ed}\n";
}
 
print "Multimap profile\n";
my %multimap_tally = ();
foreach my $mapped_time (sort {$a <=> $b} keys %multimappers){
   #print "$region\t$region{$region}\n";
   my $count = $multimappers{$mapped_time};
   if ($mapped_time <= 10){
      $multimap_tally{$mapped_time} += $count;
   } else {
      $multimap_tally{'42'} += $count;
   }
}
 
foreach my $mapped (sort {$a <=> $b} keys %multimap_tally){
   my $tally = $multimap_tally{$mapped};
   if ($mapped == 1){
      print "Mapped uniquely\t$tally\n";
   } elsif ($mapped == 42){
      print "Mapped to more than 10 location\t$tally\n";
   } else {
      print "Mapped to $mapped locations\t$tally\n";
   }
}
 
exit(0);
