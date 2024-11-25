#!/usr/bin/perl
use strict;
use warnings;

#Usage:Position to RsID

my $file=$ARGV[0];#SNP Œª÷√¡–±Ì CMP_cat.stat
my $rs=$ARGV[1];#snp151.txt

my $out=$ARGV[2];#.rsID

open(F,$rs);my %rs;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $pos=$a[1].":".$a[2];
	$rs{$pos}=$a[4];

}
close F;
print "Start converting\n";
open(F,$file);
open(OUT,">$out");my $num=0;my $no=0;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	if(exists $rs{$a[1]}){$num++;}else{$no++;}
	print OUT "$l\t$rs{$a[1]}\n";
}
close F;
print "Finish convert\n";
print "Rs finder:$num\nNo finder:$no\n";