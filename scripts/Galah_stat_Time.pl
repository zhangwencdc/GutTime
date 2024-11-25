#!/usr/bin/perl
use strict;
#use warnings;
use File::Basename qw(basename dirname);
#Usage:Galah½á¹û·ÖÎö

my $galah=$ARGV[0];#clusters.tsv
my $anno=$ARGV[1];#gtdbtk.bac120.summary.tsv
my $stat=$ARGV[2];#China_Bins.seqkit
my $time=$ARGV[3];#ID.txt

open(F,$time);my %people;my %time;

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split "\t",$l;
	$people{$a[0]}=$a[1];
	$time{$a[0]}=$a[2];
}
close F;


open(F,$stat);my %stat;

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split " ",$l;
	my $c=basename($a[0]);
	$stat{$c}=$a[4];#print "$c,$a[4]\n";
}
close F;

open(F,$anno);my %anno;

while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split "\t",$l;
	$anno{$a[0]}=$a[1];
}
close F;

open(FILE,$galah);my %cluster;my %cn;my %cname;
while(1){
	my $l=<FILE>;
	unless($l){last;}
	chomp $l;
	my @a=split "\t",$l;
	my $c=basename($a[0]);
	my $n=basename($a[1]);
	$cluster{$n}=$c;
	$cn{$c}++;
	if(exists $cname{$c}){$cname{$c}.=",".$n;}else{$cname{$c}=$n;}
	
}
close FILE;

open(OUT,">$ARGV[4]");
open(O,">$ARGV[4].time");
my @c=sort keys %cn;
print OUT "Cluster\tGenomeLength\tBins_Num\tTaxonomy\tDetail\n";
foreach my $c (@c) {
	my $name=substr($c,0,length($c)-4);
	my $d=basename($c);
	print OUT "$c\t$stat{$d}\t$cn{$c}\t$anno{$name}\t$cname{$c}\t";
	my $type=0;
	if($cn{$c}>=3){
		my %pv;
		my @d=split",",$cname{$c};my %tmin;my %tmax;
		foreach my $d (@d) {
			my @e=split"\\.",$d;my $e=$e[0];
			$pv{$people{$e}}++;
			my $time=$time{$e};
			if(exists $tmin{$people{$e}}){
				my $tmin=$tmin{$people{$e}};
				my @tm=split"_",$tmin;
				my @tnew=split"_",$time;
				if($tnew[0]<$tm[0]){
					$tmin{$people{$e}}=$time;
				}
				if($tnew[0]==$tm[0] && $tnew[1]<$tm[1]){
					$tmin{$people{$e}}=$time;
				}
			}else{
				$tmin{$people{$e}}=$time;
			}
			if(exists $tmax{$people{$e}}){
				my $tmax=$tmax{$people{$e}};
				my @tm=split"_",$tmax;
				my @tnew=split"_",$time;
				if($tnew[0]>$tm[0]){
					$tmax{$people{$e}}=$time;
				}
				if($tnew[0]==$tm[0] && $tnew[1]>$tm[1]){
					$tmax{$people{$e}}=$time;
				}
			}else{
				$tmax{$people{$e}}=$time;
			}
		}
		my @pv=sort keys %pv;
		foreach my $pv (@pv) {
			print OUT "$pv:$pv{$pv} $tmin{$pv} $tmax{$pv};";
			if($pv{$pv}>3){$type=1;

			my @t1=split"_",$tmin{$pv};
			my @t2=split"_",$tmax{$pv};
			my $tlength=($t2[0]-$t1[0])*12+($t2[1]-$t1[1]);

			print O "$c\t$stat{$d}\t$anno{$name}\t$pv\t$pv{$pv}\t$tmin{$pv}\t$tmax{$pv}\t$tlength\n";
			}
		}
	}
	if($type>0){print OUT "\tMultiple-timepoint";}
	print OUT "\n";
}
