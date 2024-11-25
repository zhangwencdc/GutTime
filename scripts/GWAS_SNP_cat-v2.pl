#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#Usage:GWAS SNP结果合并以及统计

my $file=$ARGV[0];#CMP_Species
my $prefix=$ARGV[1];# .FarmCPU.csv
my $rs=$ARGV[2];#snp151.txt  注意必须与比对的human参考库版本一致，这里统一用hg38版
my $out=$ARGV[3];#CMP_species.cat.v2

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

open(F,$rs);my %rs;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $pos=$a[1].":".$a[3];
	$rs{$pos}=$a[4];

}
close F;

my @file=glob "$file/*$prefix";
open(O,">$out");my %snp;my %type;my %pos;
print O "ID,Bacteria,CHROM,POS,REF,ALT,Effect,SE,Signals,Chromome,Pos,RsID\n";my $snpid=0;
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
		unless($a[7]=~/[0-9]/){next;}#仅统计有signal的snp
		unless(length($a[3]) == length($a[4])){next;}#不统计indel
		$snpid++;
		print O "$snpid,$name,$l,";
		$a[0]=~s/"//g;
		$snp{$a[0]}++;
		$type{$name}{$a[0]}++;
		
		my @b=split"-",$a[0];
		print O "$ch{$b[0]},$b[1],";$pos{$a[0]}=$ch{$b[0]}.":".$b[1];
		if(exists $rs{$pos{$a[0]}}){print O "$rs{$pos{$a[0]}}\n";}else{print O "\n"}
	}
	close F;
}
close O;
open(O2,">$out.stat");my $snpt;my $snpr;
my @snp=sort keys %snp;
my @type=sort keys %type;
foreach my $snp (@snp) {
	$snp=~s/"//g;
	print O2 "$snp\t$pos{$snp}\t$snp{$snp}\t";
	foreach my $type (@type) {
		if(exists $type{$type}{$snp}){print O2 "$type,";}
	}
	$snpt++;
	if(exists $rs{$pos{$snp}}){
	print O2 "$rs{$pos{$snp}}\n";
	$snpr++;
	}else{
	print O2 "\n";
	}
}
close O2;
print "Total SNPs:$snpt\nRsID snps:$snpr\n";