#!/usr/bin/perl
use strict;
#use warnings;
use File::Basename qw(basename dirname);
#Usage:Galah½á¹û·ÖÎö

my $galah=$ARGV[0];#clusters.tsv
my $anno=$ARGV[1];#gtdbtk.bac120.summary.tsv
my $stat=$ARGV[2];#China_Bins.seqkit

open(F,$stat);my %stat;

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split " ",$l;
	my $c=basename($a[0]);
	$stat{$c}=$a[4];#print "$c,$a[4]\n";
}
close F;

open(F,$anno);my %anno;

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split "\t",$l;
	$anno{$a[0]}=$a[1];
}
close F;

open(FILE,$galah);my %cluster;my %cn;my %cname;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split "\t",$l;
	my $c=basename($a[0]);
	my $n=basename($a[1]);
	$cluster{$n}=$c;
	$cn{$c}++;
	if(exists $cname{$c}){$cname{$c}.=",".$n;}else{$cname{$c}=$n;}
}
close FILE;

open(OUT,">$ARGV[3]");

my @c=sort keys %cn;
print OUT "Cluster\tGenomeLength\tBins_Num\tTaxonomy\tDetail\n";
foreach my $c (@c) {
	my $name=substr($c,0,length($c)-4);
	my $d=basename($c);
	print OUT "$c\t$stat{$d}\t$cn{$c}\t$anno{$name}\t$cname{$c}\n";
}
