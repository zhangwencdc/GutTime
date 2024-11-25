#!/usr/bin/perl
use strict;
#use warnings;

my $file=$ARGV[0];#Time_All_sstr_data.tsv
my $id=$ARGV[1];

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
	$timep{$a[0]}=$a[3];
}
close F;

open(OUT,">$file.people");
open(FILE,$file);my %sp;my %dp;my %species;my %compare_same;my %compare_diff;my %month;my %month_share;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
	unless(substr($a[0],0,1) eq "S"){$a[0]="S".$a[0];}
	unless(substr($a[1],0,1) eq "S"){$a[1]="S".$a[1];}
	print OUT "$l\t";
	print OUT "$people{$a[0]}\t$people{$a[1]}\t$time{$a[0]}\t$time{$a[1]}\t$timep{$a[0]}\t$timep{$a[1]}\t";$species{$a[2]}++;
	my @t1=split"_",$time{$a[0]};
	my @t2=split"_",$time{$a[1]};
	my $y1=$t1[0];my $m1=$t1[1];my $y2=$t2[0];my $m2=$t2[1];
	if($y1 == $y2 && $people{$a[0]} eq $people{$a[1]}){
		my $md=$m2-$m1;
		if($md==1 || $md==-1){
			$month{$a[2]}++;
			if($l=~/shared_strain/){$month_share{$a[2]}++;}
		}
	}
	if($l=~/shared_strain/){
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP\t";$sp{$a[2]}++;$compare_same{$a[2]}++;
				
		}else{print OUT "DiffP\t";$dp{$a[2]}++;$compare_diff{$a[2]}++;}
	}else{
		if($people{$a[0]} eq $people{$a[1]}){print OUT "SameP\t";$compare_same{$a[2]}++;
			my $size=$timep{$a[0]}-$timep{$a[1]};
			if($size==1 || $size==-1){
				if($l=~/other_strain/ && $a[3]<0.99){print OUT "ChangePoint";} #相似度小于99%，比对长度超过10kb 认为是不同菌株
			}
		}else{print OUT "DiffP\t";$compare_diff{$a[2]}++;}
	}
	print OUT "\n";
}
close FILE;
open(O,">$file.people.stat");
my @species=sort keys %species;
print O "Species,Compare within SampePeople,Share within Samepeople,Compare Between Different People,Share Between Different People,Compare within SampePeople within a month,Share within Samepeople within a month\n";
foreach my $species (@species) {
	print O "$species,$compare_same{$species},$sp{$species},$compare_diff{$species},$dp{$species},$month{$species},$month_share{$species}\n";
}