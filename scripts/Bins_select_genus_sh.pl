#!/usr/bin/perl
use strict;
use warnings;



my $file="Species.list";
my $outdir="/home/zhangwen/project/2023Time/Genome/Time/Genus";
open(F,$file);my %genus;
while(1){
	my $line=<F>;
	unless($line){last;};
	chomp $line;
	unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
	my @a=split"_",$line;
	if(exists $genus{$a[0]}){next;}
	my $genus=$a[0];$genus{$genus}++;
	system "mkdir $outdir/$genus";
	system "perl Bins_select_genus.pl gtdbtk.bac120.summary.tsv.filter $genus /home/zhangwen/project/2023Time/Genome/Time/Bins/\n";
}
close F;line){last;}
	chomp $line;
	my @a=split"\t",$line;

	my @str=split"s__",$a[1];
	my $sp=pop @str;
	$sp{$a[0]}=$sp;
}
close I;

open(F,$file);my %genus;
while(1){
	my $line=<F>;
	unless($line){last;};
	chomp $line;
	unless(substr($line,length($line)-1,1)=~/[0-9a-zA-Z]/){$line=substr($line,0,length($line)-1);}
	my @a=split"_",$line;
	if(exists $genus{$a[0]}){next;}
	my $genus=$a[0];$genus{$genus}++;
#	system "mkdir $outdir/$genus";
	system "perl Bins_select_genus.pl $info $genus /home/zhangwen/project/2023Time/Genome/Time/Bins/\n";  #挑选Bins，计算ANI
	
	#统计是否有相同菌株
	my $ANI=$outdir."/".$genus."/".$genus.".ANI";
	open(A,$ANI);my $strain=1;my %strain;my %t;my %gp;
	open(O,">$ANI.strain");
	while(1){
		my $l=<A>;
		unless($l){last;}
		chomp $l;
		my @b=split"\t",$l;
		if($b[0] eq $b[1]){next;}
		if($b[2]>=99.95){
			print O "$l\tSameStr\t";
			
		}elsif($b[2]<98.3){
			print O "$l\tDiffStr\t";
		}else{
			print O "$l\tNA\t";
		}
		$b[0] =~ /\/([A-Z0-9]+)\./;
			my $n1=$1;
			$b[1] =~ /\/([A-Z0-9]+)\./;
			my $n2=$1;
			my $tm0=$b[0];
			my $tm1=$b[1];
			$tm0=~ /\/(\S+\.fa)/;
			my $sp0=$1;
			$tm1=~ /\/(\S+\.fa)/;
			my $sp1=$1;
			print O "$sp{$sp0}\t$p{$n1}\t$time{$n1}\t$tp{$n1}\t$sp{$sp1}\t$p{$n2}\t$time{$n2}\t$tp{$n2}\t";
			if($p{$n1} eq $p{$n2}){print O "SameP\n";
			# if($b[2]>=99.95){
				if(exists $gp{$p{$n1}} ){
					if( $gp{$p{$n1}} ne $sp{$sp0}){
						print "$genus,$p{$n1},$gp{$p{$n1}},$sp{$sp0}\n";
					}
				}else{
					$gp{$p{$n1}}=$sp{$sp0};
				}
			# }
			}else{print O "DiffP\n";}

	}
	close A;
	close O;

}
close F;