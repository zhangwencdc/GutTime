#!/usr/bin/perl
use strict;
#use warnings;
#Author: Wen Zhang;2024-8
#Usage:仅保留align>10kb且来源于相同志愿者不同采样时间的样本间的比较结果

my $p=$ARGV[0];#Pairname;
my $tax=$ARGV[1];#metaphlan.taxonomy
my $f=$ARGV[2];#sstr_strain_events.tsv
my $out=$ARGV[3];#sstr_strain_events.tsv.filter

open(F,$p);
my %people;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$people{$a[2]}=$a[0];
	$people{$a[3]}=$a[0];
}
close F;

open(F,$tax);
my %tax;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $t=pop @a;my $s=pop @a;
	$tax{$t}=$s;
}
close F;

open(F,$f);
open(OUT,">$out");my %compare;my %same;my %diff;my %sp;my %same_sp;my %diff_sp;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
	unless(exists $people{$a[0]} && exists $people{$a[1]} && $people{$a[0]} eq $people{$a[1]}){next;}
	print OUT "$people{$a[0]}\t$tax{$a[2]}\t$l\n";
	$compare{$people{$a[0]}}++;$sp{$tax{$a[2]}}++;
	if($l=~/shared_strain/){$same{$people{$a[0]}}++;$same_sp{$tax{$a[2]}}++;}
	if($a[3]<0.995){$diff{$people{$a[0]}}++;$diff_sp{$tax{$a[2]}}++;}
}
close F;

my @p=keys %compare;
print OUT "People\tCompare\tSame\tDiff\n";
foreach my $p (@p) {
	print OUT "$p\t$compare{$p}\t$same{$p}\t$diff{$p}\n";
}

my @sp=sort keys %sp;
print OUT "\nSpecies\tCompare\tSame\tDiff\n";
foreach my $sp (@sp) {
	print OUT "$sp\t$sp{$sp}\t$same_sp{$sp}\t$diff_sp{$sp}\n";
}