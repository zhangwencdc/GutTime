#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#sstr_data.tsv
my $id=$ARGV[1];#id.txt
my $out=$ARGV[2];
my %people;my %time;my %timep;
open(F,$id);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];
	$time{$a[0]}=$a[2];
	$timep{$a[0]}=$a[3];
}
close F;

open(OUT,">$out");
open(FILE,$file);my %compare;my %diff;my %t;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
	unless(substr($a[0],0,1) eq "S"){$a[0]="S".$a[0];}
	unless(substr($a[1],0,1) eq "S"){$a[1]="S".$a[1];}
	unless($people{$a[0]} eq $people{$a[1]}){next;}
	
	my @t1=split"_",$time{$a[0]};
	my @t2=split"_",$time{$a[1]};
	my $y1=$t1[0];my $m1=$t1[1];my $y2=$t2[0];my $m2=$t2[1];
	my $t;
	if($y1<$y2){
		$t=($y2-$y1)*12+$m2-$m1;
	}
	if($y1==$y2){
		$t=$m2-$m1;
		if($t<0){$t=-$t;}
	}
	if($y1>$y2){
		$t=($y1-$y2)*12+$m1-$m2;
	}
	
			print OUT "$l\t$t\t";$t{$t}++;
			print OUT "$people{$a[0]}\t$people{$a[1]}\t$time{$a[0]}\t$time{$a[1]}\t$timep{$a[0]}\t$timep{$a[1]}\n";
			my $n=$a[0]."\t".$a[1];
			$compare{$n}{$t}++;
			if($a[3]<0.995){$diff{$n}{$t}++;}

	
	
}
close FILE;

my @key=keys %compare;my @t=sort keys %t;my %rate1;my %rate6;
print OUT "Sample1\tSample2\t";
open(O,">$out.v2");
foreach my $t (@t) {
	print OUT "$t M\t";
}
foreach my $n (@key) {
	print OUT "$n\t";
	foreach my $t (@t) {
		print OUT "$compare{$n}{$t}:$diff{$n}{$t}\t";
		if($t==1 && $compare{$n}{$t}>0){
			my $r=$diff{$n}{$t}/$compare{$n}{$t};
			print O "1M\t$n\t$r\n";
		}
		if($t==6 && $compare{$n}{$t}>0){
			my $r=$diff{$n}{$t}/$compare{$n}{$t};
			print O "6M\t$n\t$r\n";
		}
	}
	print OUT "\n";
}