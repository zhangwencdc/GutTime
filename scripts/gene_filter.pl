#!/usr/bin/perl
use strict;
use warnings;

#Usage:�޳�gene���Ӽ�������g__��s__���У�

open(F,"Time_Humann-cat-genefamilies-cpm.csv");
open(OUT,">Time_Humann-cat-genefamilies-cpm.csv.filter");

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	if($a[0]=~/g__/ || $a[0]=~/s__/ || $a[0]=~/unclassified/){next;}
	print OUT "$l\n";
}
close F;