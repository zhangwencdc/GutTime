#!/usr/bin/perl
use strict;
#use warnings;

#conda activate samestr


my $npy="/home/zhangwen/project/2022Time/WGS_Analysis/Metaphlan/Time/Samestr/samestr/All_species";
my $marker_db="/home/zhangwen/Data/Samestr";
#my $outdir=$ARGV[0]; #/home/zhangwen/project/2022Time/WGS_Analysis/Samestr
my $id=$ARGV[0];#ID_2023.txt


my @species=glob "$npy/*.npy";

foreach my $sp (@species) {
	my $line=substr($sp,0,length($sp)-4);
	if(-e "$line/sstr_data.tsv"){next;}
	system "samestr filter --input-files $sp --input-names $line.names.txt --marker-dir $marker_db/ --samples-min-n-hcov 5000 --species-min-samples 2 --marker-trunc-len 20 --global-pos-min-n-vcov 2 --sample-pos-min-n-vcov 5 --sample-var-min-f-vcov 0.1 --nprocs 4 --output-dir $line\n";
	#if(-e "$line/$line.npy"){
		system "samestr compare --input-files $line/*.npy --input-names $line/*.names.txt --nprocs 4 --output-dir $line/com/\n";
		system "samestr summarize --input-dir $line/com/ --mp-profiles-dir /home/zhangwen/project/2022Time/WGS_Analysis/Metaphlan/Time/Samestr/ --output-dir $line\n";
		system "samestr stats --input-files $line/*.npy --input-names $line/*.names.txt --nprocs 4 --output-dir $line/stats \n";
#	}else{
#		system "rm -rf $line\n";
#}
system "rm -rf tmp\n";
}

system "cat */sstr_data.tsv >Cat_sstr_data.tsv\n";

#Change ID
my $file="Cat_sstr_data.tsv";

my %id;my %people;my %time;
open(F,$id);
while(1){
	my $l=<F>;
	unless($l){last;}
	chomp $l;
	unless(substr($l,length($l)-1,1)=~/[0-9a-zA-Z]/){$l=substr($l,0,length($l)-1);}
	my @a=split"\t",$l;
	$id{$a[0]}=$a[3];
	$people{$a[0]}=$a[1];
	$time{$a[0]}=$a[2];
}
close F;

open(FILE,$file);my %mut;
my $line=<FILE>;
chomp $line;
my @name=split"\t",$line;
open(OUT,">$file.changeID");
while(1){
	my $line=<FILE>;
	unless($line){last;}
	chomp $line;

	my @a=split"\t",$line;
	my $num=@a;
	print OUT "$line\t$id{$a[0]}\t$people{$a[0]}\t$time{$a[0]}\t$id{$a[1]}\t$people{$a[1]}\t$time{$a[1]}\n";

}
close FILE;
close OUT;

###Stat_1month.pl
open(F,"$file.changeID"); my %share;my %not_share;my %p;my %species;
while(1){
	my $line=<F>;
	unless($line){last;}
	chomp $line;
	my @a=split"\t",$line;
	unless($a[4]>=10000){next;}  #仅统计overlap大于10kb的结果
	unless($a[7] eq "other_strain" || $a[7] eq "shared_strain"){next;}
	unless($a[9] eq $a[12]){next;} #仅统计相同个体来源的数据;
	my $t1;my $t2;
	if(substr($a[10],0,1) eq "T"){$t1=substr($a[10],1);}else{$t1=$a[10];}
	if(substr($a[13],0,1) eq "T"){$t2=substr($a[13],1);}else{$t2=$a[13];}
	my $inter=$t1-$t2;
	if($inter<0){$inter=-$inter;}
	unless($inter<=1){next;}
	if($a[7] eq "shared_strain"){$share{$a[9]}{$a[2]}++;}
	if($a[7] eq "other_strain"){$not_share{$a[9]}{$a[2]}++;}
	$p{$a[9]}++;$species{$a[2]}++;
}
close F;

my @key=sort keys %p; my %t_share;my %t_notshare; my @sp=sort keys %species;

print "Stat Doing\n";
open(OUT,">$file.changeID.stat");
print OUT "Species\tPeople\tShare\tNot_share\n";
foreach my $sp (@sp) {
	my $sum_s;my $sum_d;
	foreach my $k (@key) {
		print OUT "$sp\t$k\t$share{$k}{$sp}\t$not_share{$k}{$sp}\t";
		$sum_s+=$share{$k}{$sp};
		$sum_d+=$not_share{$k}{$sp};
		if(exists $not_share{$k}{$sp}){
			my $dis=$not_share{$k}{$sp}/($not_share{$k}{$sp}+$share{$k}{$sp});
			print OUT "$dis\n";
		}else{print OUT "\n";}
		
	}
	print OUT "$sp\tAllPeople\t$sum_s\t$sum_d\t";
	my $a=$sum_s+$sum_d;
	my $v;
	if($a>0){$v=$sum_d/$a;}
	print OUT "$v\n";
}