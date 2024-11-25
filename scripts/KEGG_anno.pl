#!/usr/bin/perl
use strict;
use warnings;

my $file=$ARGV[0];
my $out=$ARGV[1];

my $anno="kegg_hierarchy.txt";
my %anno;
open(A,$anno);
while(1){
	my $line=<A>;
	unless($line){last;}
	chomp $line;
#	$line=~tr/\"//g;
	print $line;
	my @a=split"\t",$line;
	$anno{$a[1]}=$line;
}
close A;

open(FILE,$file);
open OUT, ">$out";
my $line=<FILE>;
print OUT "$line\n";
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	unless($a[0]=~/^ko/){next;}
	if($a[0]=~/g__/ || $a[0]=~/unclassified/){
		next;
	}
	print OUT "$line\t$anno{$a[0]}\n";
}
close FILE;