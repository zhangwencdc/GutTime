#!/usr/bin/perl
use strict;
use warnings;

#Usage:查找参加两次项目的志愿者，并列出其编号
#perl Pairname.pl ProjectA.list ProjectB.list Pairname
my $f1=$ARGV[0];#ProjectA.list
my $f2=$ARGV[1];#ProjectB.list
my $out=$ARGV[2];#Pairsample

open(F,$f1);my %p;my %phone;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $p=$a[0];
	my $phone=pop @a;
	if($a[1]=~/[0-9a-zA-Z]/){
	$p{$a[0]}=$a[1];
	$phone{$phone}=$a[1];
	}
}
close F;

open(F,$f2);
open(OUT,">$out");
print OUT "Name\tPhone\tProject1\tProject2\n";
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	my $p=$a[0];
	my $phone=pop @a;
	if(exists $p{$p}){
		if($a[1]=~/[0-9a-zA-Z]/){
			print OUT "$p\t$phone\t$p{$p}\t$a[1]\n";
		}
	}elsif(exists $phone{$phone}){
		if($a[1]=~/[0-9a-zA-Z]/){
			print OUT "$p\t$phone\t$phone{$phone}\t$a[1]\n";
		}
	}
}
close F;
