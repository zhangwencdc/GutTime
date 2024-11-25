#!/usr/bin/perl
use strict;
use warnings;

my $core=$ARGV[0];#result/coregene.csv
my $anno=$ARGV[1];#prokka/Ref.gff
my $out=$ARGV[2];#Akk_coregene.anno

open(F,$core);my %cg;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split",",$l;
	$cg{$a[1]}=$l;
}
close F;
my @cg=sort keys %cg;
open(OUT,">$out");
open(FILE,$anno);
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	foreach my $cg (@cg) {
		if($line=~/$cg/){
			print OUT "$cg{$cg},$line\n";
		}
	}
}
close FILE;
close OUT;