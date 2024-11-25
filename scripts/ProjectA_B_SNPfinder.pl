#!/usr/bin/perl
use strict;
#use warnings;

#Usage：64个特异性SNP在ProjectA_B/CMP_6M样本中的存在情况

my $snp=$ARGV[0];#CMP_species.cat.v3.cut100
my $target=$ARGV[1];#ProjectA_B_GWAS.merge.vcf.filter

my %snp;

open(F,$snp);
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split",",$line;
	$a[2] =~ s/"//g;
	if($a[9]>1e-50){next;}#仅考虑signal小于1e-5的SNP
	$snp{$a[2]}=$a[9];
	#print "$a[2]\n";
}
close F;

my @snp=keys %snp;

open(OUT,">ProjectA_B_SNP.finder");
open(F,$target);my $match=0;
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	#print "$a[0],$a[1]\n";
	if($a[0]=~/CHROM/){print OUT "$line\tSignal\n";next;}
	my $i=$a[0]."-"	.$a[1];
	if(exists $snp{$i}){print OUT "$line\t$snp{$i}\n";$match++;
	my $count1 = () = $line =~ /0\/0/g;
	my $count2 = () = $line =~ /1\/1/g;
	print "$i,ref:$count1,snp:$count2\n";
	}

}
close F;

print "Match SNPs:$match\n";