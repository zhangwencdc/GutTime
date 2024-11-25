#!/usr/bin/perl
use strict;
use warnings;

#conda activate biobakery

my $file="Species.list";
my $outdir="/home/zhangwen/project/2022Time/WGS_Analysis/Strainphlan";
open(F,$file);
while(1){
	my $line=<F>;
	unless($line){last;};
	chomp $line;
	unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
	system "mkdir $outdir/$line";
	system "strainphlan -m /home/zhangwen/Data/Metaphlan/marker_db/$line.markers.fa -s /home/zhangwen/project/2022Time/WGS_Analysis/Metaphlan/Time/StrainPhlan/consensus_markers/*.pkl -o $outdir/$line/  -c $line --mutation_rates\n";
}
close F;