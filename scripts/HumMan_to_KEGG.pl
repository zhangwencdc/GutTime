#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#��HumannĬ�Ͻ��ת��ΪKEGG

#/home/zhangwen/project/2022Time/WGS_Analysis/Metaphlan/Time/HuManN-v2/KEGG

my $file=$ARGV[0];  ##�����ļ�·��  ΪHumann3������·�����������*_genefamilies.tsv *_pathabundance.tsv��*_pathcoverage.tsv����
my $out=$ARGV[1];#������·��

my @file=glob "$file/*_genefamilies.tsv.cpm";

foreach my $f (@file) {
	my $name=basename($f);
	my $out1=$out."/".$name."_kegg_genefamilies.tsv";
	system "humann_regroup_table --input $f --groups uniref90_rxn --output $out1\n";
	system "humann3 -i $out1 -o $out --pathways-database /home/zhangwen/Data/Metaphlan/HuManN/mapping_v201901/keggc\n";
}