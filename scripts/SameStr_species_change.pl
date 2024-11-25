#!/usr/bin/perl
use strict;
#use warnings;

#Usage:统计249个Species（在三个以上样本中发现的菌种）是否在本研究中发生了菌株水平上的变化

my $f=$ARGV[0];#Species.list
my $file=$ARGV[1];#sstr_data.tsv.people
my $out=$ARGV[2];#sstr_data.tsv.people.change
my %species;
open(F,$f);
while(1){
	my $l=<F>;
	unless($l){last;}
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	$species{$l}++;
}
close F;

open(FILE,$file);my %change;my %compare;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	my @a=split"\t",$l;
	unless($l=~/SameP/){next;}
	unless(exists $species{$a[2]}){next;}
	$compare{$a[2]}++;
	unless($l=~/other_strain/&& $a[3]<0.995){next;}
	#print "$l\n";
	$change{$a[8]}{$a[2]}++
}
close FILE;

my @p=sort keys %change;
my @s=keys %species;
open(OUT,">$out");
print OUT "Species\t";
foreach my $p (@p) {
	print OUT "$p\t";
}
print OUT "Sum Compare\n";
foreach my $s (@s) {
	print OUT "$s\t";
	foreach my $p (@p) {
	print OUT "$change{$p}{$s}\t";
	}
	print OUT "$compare{$s}\n";
}