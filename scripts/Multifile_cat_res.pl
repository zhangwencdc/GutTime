#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);

#Usage:合并bowtie比对的多文件结果 至一个文件
#Example:perl Multifile_cat_res.pl .res.bam.sort.coverage Time_res.csv

my $file=$ARGV[0];# .res.bam.sort.coverage

my @file=glob "./*$file";

my %gene;my %anno;my %sample;
foreach my $file (@file) {
	my $name=basename($file);
	open(F,$file);
	$sample{$name}=0;
	while(1){
		my $line=<F>;
		unless($line){last;}
		chomp $line;
		if($line=~/^Gene ID/){next;}
		my @a=split"\t",$line;
		#my $align=$a[3]/$a[1];
		unless($a[5]>=90){next;}
		$sample{$name}++;
		$gene{$a[0]}{$name}=$a[5];
		#$anno{$a[0]}=$a[6];
	}
	close F;

}
my $out=$ARGV[1];
open(OUT,">$out");
my @sample=sort keys %sample;

print OUT "ID,";
print "Sample,Match Num\n";
foreach my $sample (@sample) {
	print OUT "$sample,";
	print "$sample,$sample{$sample}\n";
}
print OUT "\n";

my @gene=sort keys %gene;
foreach my $gene (@gene) {
	print OUT "$gene,";
	foreach my $sample (@sample) {
		if(exists $gene{$gene}{$sample}){
				print OUT "$gene{$gene}{$sample},";
		}else{
			print OUT "0,";
		}
	}
	#my $anno=$anno{$gene};
	print OUT "\n";
}
