#!/usr/bin/perl
use strict;
use warnings;

my $file=$ARGV[0];#China_sstr_data_samestrain

my %sp;my %strain; my %type;
open(F,$file);
my %com;
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	my @a=split"\t",$l;
	unless($a[4]>=10000){next;}
	unless($a[3]>=0.9999){next;}
	unless($l=~/shared_strain/){next;}
	my $t1=$a[0]."_".$a[2];
	my $t2=$a[1]."_".$a[2];
	my $pc=$a[0]."_".$a[1];
	$com{$pc}++;
	if(exists $strain{$a[0]}{$a[2]}){
		if(exists $strain{$a[1]}{$a[2]}){
			if($strain{$a[0]}{$a[2]}<$strain{$a[1]}{$a[2]}){
				$strain{$a[1]}{$a[2]}=$strain{$a[0]}{$a[2]};
				my @type=keys %type;
				foreach my $type (@type) {
					$type{$t2}=$type{$t1};
				}
			}else{
				$strain{$a[0]}{$a[2]}=$strain{$a[1]}{$a[2]};
				my @type=keys %type;
				foreach my $type (@type) {
					$type{$t1}=$type{$t2};
				}
			}
		}else{
		$strain{$a[1]}{$a[2]}=$strain{$a[0]}{$a[2]};
		$type{$t2}=$type{$t1};
		}
	}elsif(exists $strain{$a[1]}{$a[2]}){
		if(exists $strain{$a[0]}{$a[2]}){
			if($strain{$a[0]}{$a[2]}<$strain{$a[1]}{$a[2]}){
				$strain{$a[1]}{$a[2]}=$strain{$a[0]}{$a[2]};
				my @type=keys %type;
				foreach my $type (@type) {
					$type{$t2}=$type{$t1};
				}
			}else{
				$strain{$a[0]}{$a[2]}=$strain{$a[1]}{$a[2]};
				my @type=keys %type;
				foreach my $type (@type) {
					$type{$t1}=$type{$t2};
				}
			}
		
		}else{
			$strain{$a[0]}{$a[2]}=$strain{$a[1]}{$a[2]};
			$type{$t1}=$type{$t2};
		}
	}else{
		$sp{$a[2]}++;
		$strain{$a[1]}{$a[2]}=$sp{$a[2]};
		$strain{$a[0]}{$a[2]}=$sp{$a[2]};
		$type{$t1}=$a[2]."_".$strain{$a[0]}{$a[2]};
		$type{$t2}=$a[2]."_".$strain{$a[0]}{$a[2]};
	}
}
close F;

my @type=keys %type;my %count;
foreach my $type (@type) {
	print "$type\t$type{$type}\n";
	$count{$type{$type}}++;
}

print "\n 多个样本中都存在的stain\n";
my @count=sort  {$count{$b}<=>$count{$a}} keys %count ;
foreach my $count (@count) {
	if($count{$count}>=5){
		print "$count\t$count{$count}\n";
	}
}


print "People共享菌种数\n";
my @com=sort  {$com{$b}<=>$com{$a}} keys %com ;

foreach my $com (@com) {
	print "$com\t$com{$com}\n";
}