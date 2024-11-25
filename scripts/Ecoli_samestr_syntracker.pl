#!/usr/bin/perl
use strict;
use warnings;

#Usage：结合Syntracker和Samestr skani结果，查找Ecoli相同株

my $f1=$ARGV[0];#Syntracker 结果 avg_synteny_scores_all_regions.csv

my $f2=$ARGV[1];#Samestr_China_sstr_data.tsv.Ecoli

my $f3=$ARGV[2];#Escherichina_Bins.ANI
my $out=$ARGV[3];#Ecoli_samestr_syntracker.txt

open(F,$f3); my %ani;
while(1){
my $l=<F>;
unless($l){last;}
chomp $l;

my @a=split"\t",$l;
$a[0] =~ /(.*?)\./;my $n1=$1;
$a[1] =~ /(.*?)\./;my $n2=$1;
$ani{$n1}{$n2}=$a[2];
}

close F;

open(F,$f1); my %syntracker;
while(1){
my $l=<F>;
unless($l){last;}
chomp $l;
$l =~ s/"//g; 
my @a=split",",$l;
$syntracker{$a[1]}{$a[2]}=$a[3];
}

close F;


open(F,$f2);
open(OUT,">$out");print OUT "Sample1\tSample2\tSamestr\tSyntracker\tGenome ANI\n";
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>10000){next;}
	my @b=split"_",$a[0];
	my $n1=$b[0];
	my @c=split"_",$a[1];
	my $n2=$c[0];
	if(exists $syntracker{$n1}{$n2}){
		print OUT "$n1\t$n2\t$a[3]\t$syntracker{$n1}{$n2}\t";
	}elsif(exists $syntracker{$n2}{$n1}){
		print OUT "$n2\t$n1\t$a[3]\t$syntracker{$n2}{$n1}\t";
	}else{
		print OUT "$n2\t$n1\t$a[3]\t\t";
	}
	if(exists $ani{$n1}{$n2}){
	print OUT "$ani{$n1}{$n2}\n";
	}elsif(exists $ani{$n2}{$n1}){
	print OUT "$ani{$n2}{$n1}\n";
	}else{print OUT "\n";}

}
close F;
