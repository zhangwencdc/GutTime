#!/usr/bin/perl
use strict;
#use warnings;
#Usagea:计算样本两两间的bray距离，仅保留间隔1月的结果，并标记该月是否有抗生素用药史
my $f=$ARGV[0];#otutab.txt
my $meta=$ARGV[1];#metadata.tsv
#my $res=$ARGV[2];#抗生素用药史；
my $out=$ARGV[2];#菌种流失率.6M
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



open(F,$f);
my $l=<F>;chomp $l;
my @sample=split"\t",$l;my $n=@sample;my %dis;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	foreach  (1..($n-1)) {
		#if(exists $dis{$sample[$_]}{$[0]}){next;}
		$dis{$a[0]}{$sample[$_]}=$a[$_];
	}
}
close F;
shift @sample;my @species=keys %dis;
open(OUT,">$out");
print OUT "Sample1\tPeople1\tTime1\tSpecies1\tSample2\tPeople2\tTime2\tSpecies2\tShare_Species\n";
foreach my $s1 (@sample) {
	foreach my $s2 (@sample) {
		if($s1 eq $s2){next;}
		unless($people{$s1} eq $people{$s2}){next;}
		my $t1=$time{$s1};my $t2=$time{$s2};
		my @t1=split"_",$t1;my @t2=split"_",$t2;
		my $t;
		if($t1[0] eq $t2[0]){
			if(($t2[1]-$t1[1])==1){$t=1;}
		}elsif($t1[0] +1== $t2[0]){
			if($t1[1]==12 && $t2[1]==1){$t=1;}
		}
		unless($t==1){next;}
		my $num1;my $num2;my $num_s;
		foreach my $s (@species) {
			if($dis{$s}{$s1}>=0.01){
				$num1++;
			}
			if($dis{$s}{$s2}>=0.01){
				$num2++;
			}
			if($dis{$s}{$s1}>=0.01 && $dis{$s}{$s2}>=0.01){
				$num_s++;
			}
		}
		print OUT "$s1\t$people{$s1}\t$time{$s1}\t$num1\t$s2\t$people{$s2}\t$time{$s2}\t$num2\t$num_s\n";
	}
}