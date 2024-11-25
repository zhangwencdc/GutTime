#!/usr/bin/perl
use strict;
use warnings;

#Usage:比较Time和ProjectA_B两个队列里的菌变迁率

my $f1=$ARGV[0];#Time_rate.txt
my $f2=$ARGV[1];#ProjectA_B_6M_rate.txt
my $out=$ARGV[2];#Compare_ProjectA_B_Time

open(F,$f1);
my %rate;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$rate{$a[0]}=$a[4];
}
close F;


open(F,$f2);
open(OUT,">$out");
print OUT "Species\tTime_6Mrate\tProjectA_B_6Mrate\n";
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	$a[0]=substr($a[0],3);#剔除开头的s__字符
	if(exists $rate{$a[0]}){print OUT "$a[0]\t$a[4]\t$rate{$a[0]}\n";}
}
close F;