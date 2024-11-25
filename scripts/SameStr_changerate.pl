#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#sstr_data.tsv.people
my $out=$ARGV[1];#sstr_rate
open(FILE,$file);my %change1;my %change2;my %compare1;my %compare2;my %s; my %change3;my %compare3;

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
	if($time==1){
		
		$compare1{$a[2]}++;
		unless($l=~/other_strain/ && $a[3]<0.995){next;}
		#print "$l\n";
		$change1{$a[2]}++;
	}elsif($time==6){
		$compare2{$a[2]}++;
		unless($l=~/other_strain/ && $a[3]<0.995){next;}
		#print "$l\n";
		$change2{$a[2]}++;
	
	}elsif($time==12){
		if($l=~/Akkermansia_muciniphila/ ){print "$l\t$time\n";}
		$compare3{$a[2]}++;
		unless($l=~/other_strain/ && $a[3]<0.995){next;}
		#print "$l\n";
		$change3{$a[2]}++;
	
	}
}
close FILE;

open(OUT,">$out");
print OUT "Species\tChange 1month\tCompare 1month\tChange 6month\tCompare 6month\tChange 1Year\tCompare 1Year\n";
my @s=sort keys %s;
foreach my $s (@s) {
	print OUT "$s\t$change1{$s}\t$compare1{$s}\t$change2{$s}\t$compare2{$s}\t$change3{$s}\t$compare3{$s}\n";
}