#!/usr/bin/perl

use strict;
use warnings;

my $usage = "Usage: $0 <fasta.fa>\n";
my $infile = shift or die $usage;

my $a = 0;
my $c = 0;
my $g = 0;
my $t = 0;
my $other = 0;

open (IN, '<', $infile) or die "Can't open $infile: $!\n";
while (<IN>){
   chomp;
   next if (/^$/);
   next if (/^>/);

	my $seq = $_;

   for (my $i = 0; $i < length($seq); ++$i){
      my $base = substr($seq,$i,1);
      if ($base =~ /a/i){
         ++$a;
      } elsif ($base =~ /c/i){
         ++$c;
      } elsif ($base =~ /g/i){
         ++$g;
      } elsif ($base =~ /t/i){
         ++$t;
      } else {
         ++$other;
      }
   }
}
close(IN);

print "A: $a\n";
print "C: $c\n";
print "G: $g\n";
print "T: $t\n";
print "Other: $other\n";

exit(0);

