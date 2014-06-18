#!/bin/env perl
 
use strict;
use warnings;

my $usage = "Usage: $0 <number> <size>\n";
my $number = shift or die $usage;
my $size = shift or die $usage;

my $outfile = 'my_random_' . $size . '_' . $number . '.fa';

open(OUT,'>',$outfile) || die "Could not open $outfile: $!\n";
for (1 .. $number){
   my $seq = ran_seq($size);
   print OUT ">${_}_$size\n$seq\n";
}
close(OUT);
 
sub ran_seq {
   my ($length) = @_;
   my $seq = '';
   my @nuc = qw/A C G T/;
   for (1 .. $length){
      my $rand = int(rand(scalar(@nuc)));
      $seq .= $nuc[$rand];
   }
   return($seq);
}
 
exit(0);
