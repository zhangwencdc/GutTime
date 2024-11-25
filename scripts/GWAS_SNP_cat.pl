#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#Usage:GWAS SNP结果合并以及统计

my $file=$ARGV[0];#CMP_Genus
my $prefix=$ARGV[1];# .FarmCPU_signals.csv
my $out=$ARGV[2];#*.cat

open(CH,"Homo_chromosome.name");
my %ch;
while(1){
	my $l=<CH>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$ch{$a[0]}=$a[1];
}
close CH;

my @file=glob "$file/*$prefix";
open(O,">$out");my %snp;my %type;my %pos;
print O "SNP,CHROM,POS,REF,ALT,Effect,SE,Signals,Chromome,Pos\n";
foreach my $f (@file) {
	print "$f\n";
	my $name=basename($f);
	open(F,$f);
	my $l=<F>;chomp $l;#print O "$name,$l\n";
	while(1){
		my $l=<F>;
		unless($l){last;}
		chomp $l;
		my @a=split",",$l;
		print O "$name,$l,";
		$a[0]=~s/"//g;
		$snp{$a[0]}++;
		$type{$name}{$a[0]}++;
		
		my @b=split"-",$a[0];
		print O "$ch{$b[0]},$b[1]\n";$pos{$a[0]}=$ch{$b[0]}.":".$b[1];
	}
	close F;
}
close O;
open(O2,">$out.stat");
my @snp=sort keys %snp;
my @type=sort keys %type;
foreach my $snp (@snp) {
	$snp=~s/"//g;
	print O2 "$snp\t$pos{$snp}\t$snp{$snp}\t";
	foreach my $type (@type) {
		if(exists $type{$type}{$snp}){print O2 "$type,";}
	}
	print O2 "\n";
}
close O2;