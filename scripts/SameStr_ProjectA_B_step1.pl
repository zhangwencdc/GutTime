#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#conda activate biobakery


my @file=glob "/bk-machine/zhangwen/Fecal/Fecal_Data//*_1.clean.fq";

foreach my $file (@file) {
        my $f=basename($file);
        $f=substr($f,0,length($f)-11);
        my $f2=substr($file,0,length($file)-11)."_2.clean.fq";
        if( -e "$f.sam"){next;}
        system " metaphlan $file,$f2 --bowtie2db /data/zhangwen/Data/metaphlan_databases/ --input_type fastq --nproc 30 --legacy-output -t rel_ab --bowtie2out $f.bowtieout --samout $f.sam.bz2 -o $f.profile.txt \n";
}