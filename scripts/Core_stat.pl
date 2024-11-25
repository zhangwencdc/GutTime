#!/usr/bin/perl
use strict;
#use warnings;

#Usage：Core Genus or Core Species

my $file=$ARGV[0];#Metaphlan合并结果
my $id=$ARGV[1];#ID.txt
my $cutoff=$ARGV[2];#default 0.8

my %people;my %time;my %timep;my %p;
open(F,$id);my $l=<F>;chomp $l;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];$p{$a[1]}++;
	$time{$a[0]}=$a[2];
	$timep{$a[0]}=$a[3];
}
close F;

open(FILE,$file);
my $l=<FILE>;chomp $l;#剔除第一行 #mpa_v30_CHOCOPhlAn_201901

 $l=<FILE>;chomp $l; my @name=split"\t",$l; #读取表头
my $n=@name; my %pn;my $sum=$n-2;
foreach  (2..($n-1)) {
	my $name=$name[$_];
	if($name[$_]=~/.profile/){
		$name=substr($name,0,length($name)-8);
		if(exists $people{$name}){$pn{$people{$name}}++;}
	}
}
my @people=sort keys %p;

open(O,">$file.prevalence");
open(O2,">$file.prevalence.v2");
open(OUT,">$file.prevalence.core");
print O "People\tSum\t";print OUT "People\tSum\t";print O2 "People\tSum\t";
open(GO,">$file.people.prevalence");

foreach my $p(@people) {
	print O "$p\t";print OUT "$p\t";print O2 "$p\t";
}
print O "\n";print O2 "\n";print OUT "\n";my $pn=@people;my %o;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $target=$a[0];
	unless($target=~/g__/){next;}#仅统计genus or species 水平的结果
	my $asum;my %sp;
	foreach(2..($n-1)){
		if($a[$_]>0){
			my $name=$name[$_];
			$name=substr($name,0,length($name)-8);
			my $sp=$people{$name};
			$asum++;$sp{$sp}++;
		}
	}
	my $asump=$asum/$sum;
	print O "$target\t$asump\t";print O2 "$target\t$asum\t";
	print OUT "$target\t";my $t1;my $t2;
	if($asump>=$cutoff){print OUT "Core\t";$t1++;}else{print OUT "NA\t";}
	foreach my $people (@people) {
		my $pp=$sp{$people}/$pn{$people};
		print O "$pp\t";print O2 "$sp{$people}\t";$o{$target}{$people}=$sp{$people};
		if($pp>=$cutoff){print OUT "Core\t";$t2++;}else{print OUT "NA\t";}
	}
	print O "\n";print O2 "\n";
	if($t1>0){
		if($t2>=$pn){print OUT "Core_in_All";}elsif($t2>0){print OUT "Core_in_People";}
	}else{
		if($t2>0){print OUT "Core_in_People";}
	}
	print OUT "\n";

}
close FILE;

my @target=keys %o;my %stat;my $max=0;
foreach my $target (@target) {
	foreach my $p (@people) {
		if(exists $o{$target}{$p}){
			print GO "$target,$p,$o{$target}{$p}\n";
			$stat{$p}{$o{$target}{$p}}++;
			if($o{$target}{$p}>$max){$max=$o{$target}{$p};}
		}
	}
}
print GO "ID,";
foreach my $p (@people) {
	print GO "$p,";
}
print GO "\n";
foreach  (1..$max) {
	print GO "$_,";
	foreach my $p (@people) {
		print GO "$stat{$p}{$_},";
	}
	print GO "\n";
}