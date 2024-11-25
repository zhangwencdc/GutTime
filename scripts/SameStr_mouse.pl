#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#Mouse_sstr_strain_events.tsv
my $id=$ARGV[1];#Mouse_ID.txt
my $tax=$ARGV[2];#Mouse_metaphlan.taxonomy

my %people;my %time;my %timep;
open(F,$id);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$people{$a[0]}=$a[1];
	$time{$a[0]}=$a[2];
	#$timep{$a[0]}=$a[3];
	#print "$l\n$a[1]\t$a[2]\n";
}
close F;

my %tax;
open(F,$tax);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	my $id=pop @a;my $sp=pop @a;$tax{$id}=$sp;
}
close F;

open(OUT,">$file.mouse");
open(FILE,$file);my %sp;my %dp;my %species;my %compare_same;my %compare_diff;my %month;my %month_share;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
#	unless(substr($a[0],0,1) eq "S"){$a[0]="S".$a[0];}
#	unless(substr($a[1],0,1) eq "S"){$a[1]="S".$a[1];}
	print OUT "$l\t";my $species=$tax{$a[2]};
	print OUT "$species,$people{$a[0]}\t$people{$a[1]}\t$time{$a[0]}\t$time{$a[1]}\t";$species{$species}++;
	
	if($l=~/shared_strain/){
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP\t";$sp{$species}++;$compare_same{$species}++;
				
		}else{print OUT "DiffP\t";$dp{$species}++;$compare_diff{$species}++;}
	}else{
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP\t";$compare_same{$species}++;
			
		}else{print OUT "DiffP\t";$compare_diff{$species}++;}
	}

	if($people{$a[0]} eq $people{$a[1]} && $a[3]<0.995){print "$l\t$species\t$people{$a[0]}\t$people{$a[1]}\t$time{$a[0]}\t$time{$a[1]}\n";}
	print OUT "\n";
}
close FILE;
open(O,">$file.mouse.stat");
my @species=sort keys %species;
print O "Species,Compare within SampePeople,Share within Samepeople,Compare Between Different People,Share Between Different People,Compare within SampePeople within a month,Share within Samepeople within a month\n";
foreach my $species (@species) {
	print O "$species,$compare_same{$species},$sp{$species},$compare_diff{$species},$dp{$species},$month{$species},$month_share{$species}\n";
}