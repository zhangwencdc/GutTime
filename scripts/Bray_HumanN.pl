#!/usr/bin/perl
use strict;
#use warnings;
#Usagea:基于humanN预测的基因结果，计算样本两两间的bray距离，仅保留来自于相同People间隔1月的结果
my $f=$ARGV[0];#bray.csv
my $out=$ARGV[1];#bray_humanN.list


open(F,$f);
my $l=<F>;chomp $l;
my @sample=split",",$l;my $n=@sample;my %dis;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split",",$l;
	foreach  (1..($n-1)) {
		#if(exists $dis{$sample[$_]}{$[0]}){next;}
		$dis{$a[0]}{$sample[$_]}=$a[$_];
	}
}
close F;
shift @sample;
open(OUT,">$out");
foreach my $s1 (@sample) {
	foreach my $s2 (@sample) {
		if($s1 eq $s2){next;}
		my @p1=split"_",$s1;
		my @p2=split"_",$s2;
		unless($p1[0] eq $p2[0]){next;}
		my $t;
		if($p1[1] eq $p2[1]){
			if(($p2[2]-$p1[2])==1){$t=1;}
		}elsif($p1[1] < $p2[1]){
			if($p1[2]==12 && $p2[2]==1){$t=1;}
		}
		unless($t==1){next;}
		
		print OUT "$s1,$s2,$dis{$s1}{$s2}\n";
	}
}