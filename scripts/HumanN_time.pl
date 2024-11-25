#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#conda activate biobakery
#/home/zhangwen/project/2023Time/HumannN

my $file=$ARGV[0];#ID.txt
open(F,$file);
my $l=<F>;chomp $l;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	unless(substr($a[1],length($a[1])-1,1)=~/[0-9a-zA-Z]/){$a[1]=substr($a[1],0,length($a[1])-1);}
	my $test=$a[1]."/".$a[1]."_genefamilies.tsv";
	if(-e $test){next;}
	 print "seqkit sample -n 2000000 $a[2] >$a[1].fq\n";
	 print "humann --input $a[1].fq --output $a[1] --output-basename $a[1]\n";
	 print "rm -rf $a[1].fq\n";
}
close F;

close F;