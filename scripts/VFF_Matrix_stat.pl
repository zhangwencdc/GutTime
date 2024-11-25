#!/usr/bin/perl
use strict;
use warnings;

#Matrix Analysis

my $file=$ARGV[0];#VFF_matrix.txt
open(F,$file);
my $line=<F>;
my @id=split"\t",$line;
my %p;my %n;
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	my $name=$a[0];
	my $people=$a[1];$n{$people}++;
	my $num=@a;
	foreach  (4..($num-1)) {
		if($a[$_]>0){$p{$id[$_]}{$people}++;}
	}
}

my @people=sort keys %n;
my @k=sort keys %p;
open(OUT,">$file.people");
print OUT ",";
foreach my $p (@people) {
	print OUT "$p,";
}
print OUT "\n";
foreach my $id (@k) {
	print OUT "$id,";
	foreach my $p (@people) {
		my $v=$p{$id}{$p}/$n{$p};
		print OUT "$v,";
	}
	print OUT "\n";
}
close OUT;