#!/usr/bin/perl
use strict;
#use warnings;
#Usagea:计算样本两两间的bray距离，仅保留间隔1月的结果，并标记该月是否有抗生素用药史
my $f=$ARGV[0];#bray.csv
my $meta=$ARGV[1];#metadata.tsv
my $res=$ARGV[2];#抗生素用药史；
my $out=$ARGV[3];#bray_res.list
open(F,$meta);
my %people;my %time;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];
	$time{$a[0]}=$a[2];
}
close F;


open(F,$res);
my %res;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	$res{$l}++;
}
close F;

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
		unless($people{$s1} eq $people{$s2}){next;}
		my $t1=$time{$s1};my $t2=$time{$s2};
		my @t1=split"_",$t1;my @t2=split"_",$t2;
		my $t;
		if($t1[0] eq $t2[0]){
			if(($t2[1]-$t1[1])==1){$t=1;}
		}elsif($t1[0] < $t2[0]){
			if($t1[1]==12 && $t2[1]==1){$t=1;}
		}
		unless($t==1){next;}
		my $res=$people{$s2}."_".$time{$s2};
		print OUT "$s1,$people{$s1},$time{$s1},$s2,$people{$s2},$time{$s2},$dis{$s1}{$s2},$res{$res}\n";
	}
}