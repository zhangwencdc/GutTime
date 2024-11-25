#!/usr/bin/perl
use strict;
#use warnings;

#Usage:统计gtdbtk.bac120.summary.tsv

my $file=$ARGV[0];#gtdbtk.bac120.summary.tsv
my $id=$ARGV[1];#ID.txt
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


open(O,">$file.stat");open(OUT,">$file.filter");
open(F,$file);my $l=<F>;chomp $l;print OUT "$l\n";
my %genus;my %species;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	if($a[1]=~/Unclassified/){next;}  ##剔除未分类Bins
	unless($a[16]>=70){next;} #剔除aligh 小于70%的Bins
	my $b=pop @a;
	unless($b eq "N/A"){next;}
	my @b = split(/\./, $a[0]);
	unless(exists $people{$b[0]}){
	print "Error $a[0]\n";next;
	}
	print O "$people{$b[0]}\t$time{$b[0]}\t$timep{$b[0]}\t$a[0]\t$a[1]\n";print OUT "$l\n";
	my @c = split(/\;/, $a[1]);
	my $sp=pop @c;my $ge=pop @c;
	$genus{$ge}{$people{$b[0]}}++;
	$species{$sp}{$people{$b[0]}}++;
}
close F;

open(S,">$file.species.stat");
my @people=sort keys %p;
print S "Species\t";
foreach my $p (@people) {
	print S "$p\t";
}
print S "\n";
my @sp=sort keys %species;
foreach my $sp (@sp) {
	print S "$sp\t";my $ns;
	foreach my $p (@people) {
		print S "$species{$sp}{$p}\t";$ns+=$species{$sp}{$p};
	}
	print S "$ns\n";
}
close S;

open(S,">$file.genus.stat");

print S "Genus\t";
foreach my $p (@people) {
	print S "$p\t";
}
print S "\n";
my @sp=sort keys %genus;
foreach my $sp (@sp) {
	print S "$sp\t";my $ns;
	foreach my $p (@people) {
		print S "$genus{$sp}{$p}\t";$ns+=$genus{$sp}{$p};
	}
	print S "$ns\n";
}
close S;