#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#sstr_data.tsv.people
my $out=$ARGV[1];#sstr_duration_time
open(FILE,$file);
my %time;my %num;my %p;my %s;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	my @a=split"\t",$l;
	unless($l=~/SameP/){next;}
	my $t1=$a[10];
	my $t2=$a[11];
	my @t1=split"_",$t1;
	my @t2=split"_",$t2;
	my $y1=$t1[0];
	my $y2=$t2[0];
	my $m1=$t1[1];
	my $m2=$t2[1];
	$s{$a[2]}++;
	my $time=0;
	if($y1<$y2){
		$time=12*($y2-$y1)+($m2-$m1);
	}else{
		$time=12*($y1-$y2)+($m1-$m2);
	}
	if($time<0){$time=0-$time;}
	unless($l=~/shared_strain/ && $l=~/SameP/){next;}
	if($time>$time{$a[2]}{$a[8]}){$time{$a[2]}{$a[8]}=$time;}
	$num{$a[2]}{$a[8]}++;
	$p{$a[8]}++;
	
}
close FILE;

open(OUT,">$out");
my @p=sort keys %p;
my @s=sort keys %s;
print OUT "People\t";
foreach my $p (@p) {
	print OUT "$p\t";
}
print OUT "\n";
foreach my $s (@s) {
	print OUT "$s\t";my $sum;my $num;
	foreach my $p (@p) {
		if($num{$s}{$p}>=5){
			#my $avg=$time{$s}{$p}/$num{$s}{$p};
			$sum+=$time{$s}{$p};
			$num++;
			print OUT "$time{$s}{$p}\t";
		}else{
			print OUT "\t";
		}
	}
	if($num>0){
		my $avg=$sum/$num;
		print OUT "$avg\n";
	}else{
		print OUT "\n";
	}
}