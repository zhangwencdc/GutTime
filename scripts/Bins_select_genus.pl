#!/usr/bin/perl
use strict;
use warnings;
#Usage:基于gtdbtk预测结果提取指定菌种的基因组

my $file=$ARGV[0];#gtdbtk.bac120.summary.tsv gtdbtk预测结果

my $target=$ARGV[1];#Prevotella

#my @a=split"_",$target;
my $species="g__".$target;

my $dir=$ARGV[2];##Bins存放路径


my %num;
open(F,$file);
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	my $b=pop @a;
	unless($b eq "N/A"){next;}
	unless($a[1]=~/$species/){next;}

	$num{$a[0]}=$a[1];
}
close F;

my @sample=sort keys %num;
open(O,">$target/$target.list");
foreach my $sp (@sample) {
	my $b=$num{$sp};
	#my @b=split",",$b;
	my $out=$target."/".$sp;
	print O "$out\n";
		if($b=~/[0-9a-zA-Z]/){
			my $file=$dir."/".$sp.".fna";
			#system "cp $file $out\n";
		}
	
}

#system "fastANI --ql $target/$target.list --rl $target/$target.list -o $target/$target.ANI\n";

system "ANIclustermap -i $target/ -o $target/ --annotation";