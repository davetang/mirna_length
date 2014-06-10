#!/bin/env perl

use strict;
use warnings;

my $usage = "Usage: $0 <fasta.fa>\n";
my $infile = shift or die $usage;

my $seq = '';
open (IN, '<', $infile) or die "Can't open $infile: $!\n";
while (<IN>){
   chomp;
   next if (/^$/);
   next if (/^>/);
	$seq .= $_;

}
close(IN);

warn("Finished storing $infile\n");

my $a = 0; my $c = 0; my $g = 0; my $t = 0; my $other = 0;

for(my $i=0; $i<length($seq); ++$i){
	my $base = substr($seq,$i,1);
	if ($base =~ /a/i){
		++$a;
	}
	if ($base =~ /c/i){
		++$c;
	}
	if ($base =~ /g/i){
		++$g;
	}
	if ($base =~ /t/i){
		++$t;
	}
	if ($base !~ /[acgt]/i){
		++$other;
	}
}

print "A: $a\nC: $c\nG: $g\nT: $t\nOther: $other\n";

exit(0);

