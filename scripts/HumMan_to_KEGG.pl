#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw(basename dirname);
#将Humann默认结果转化为KEGG

#/home/zhangwen/project/2022Time/WGS_Analysis/Metaphlan/Time/HuManN-v2/KEGG

my $file=$ARGV[0];  ##输入文件路径  为Humann3结果输出路径，结果按照*_genefamilies.tsv *_pathabundance.tsv和*_pathcoverage.tsv保存
my $out=$ARGV[1];#结果输出路径

my @file=glob "$file/*_genefamilies.tsv.cpm";

foreach my $f (@file) {
	my $name=basename($f);
	my $out1=$out."/".$name."_kegg_genefamilies.tsv";
	system "humann_regroup_table --input $f --groups uniref90_rxn --output $out1\n";
	system "humann3 -i $out1 -o $out --pathways-database /home/zhangwen/Data/Metaphlan/HuManN/mapping_v201901/keggc\n";
}