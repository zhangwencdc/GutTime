#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#conda activate biobakery


my @file=glob "/home/zhangwen/project/2022Time/WGS_Data/*_R1.fq";

foreach my $file (@file) {
        my $f=basename($file);
        $f=substr($f,0,length($f)-6);
        my $f2=substr($file,0,length($file)-6)."_R2.fq";
        if( -e "$f.sam"){next;}
        system " metaphlan $file,$f2 --bowtie2db /home/zhangwen/Data/Metaphlan/ --input_type fastq --nproc 30 --legacy-output -t rel_ab --bowtie2out $f.bowtieout --samout $f.sam.bz2 -o $f.profile.txt \n";
}