#!/usr/bin/perl

=head1 Name

=head1 Description
#全新计算 或者 基于已有运算基础增加一株或几株新的基因组，重新构树
=head1 

  Author: zhangwen, zhangwen@icdc.cn
  Version: 1.0, Date: 2022-01-05

=head1 Usage:
全新 analysis有cds:
perl  core_pipeline-V2021-addnewgenome.pl -Ref_s Ref.seq -cds coregene.list  A.fa B.fa
增加新的genome：
perl  core_pipeline-V2021-addnewgenome.pl -Ref_s Ref.seq -cds coregene.list  A.fa B.fa -old all_snp.xls
=head1 Example
全新 analysis无CDS:
perl  core_pipeline-V2021-addnewgenome.pl -Ref_s Ref.seq   A.fa B.fa
增加新的genome：
perl  core_pipeline-V2021-addnewgenome.pl -Ref_s Ref.seq -cds coregene.list  A.fa B.fa -old all_snp.xls
减少程序依赖，目前仅需要dnadiff\prokka\blat\fasttree
=cut
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use Pod::Text;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET=1;

use lib '/home/zhangwen/miniconda3//pkgs/perl-bioperl-1.6.924-4/lib/perl5/site_perl/5.22.0/';
##get options from command line into variables and set default values
###-----------------Time record-----------------------------###
my $Time_Start = sub_format_datetime(localtime(time())); #......
my $Data_Vision = substr($Time_Start,0,10);
my ($list,$data_path,$outdir,$mark,$Ref,$Ref_s,$dir,$identity,$coverage,$type,$HELP,$NUMREPS,$BURNIN,$sum_n,$STEP,$coregene,$old);
my $dir;
GetOptions(
		"O:s"=>\$outdir,    
		"R|Ref_seq:s"=>\$Ref_s,
		"G|cds:s"=>\$coregene,
		"D|Dir|dir:s"=>\$dir,
		"old:s"=>\$old,
		"help"=>\$HELP
);
die `pod2text $0` if ($HELP  || !defined $Ref_s  );
if(!defined $outdir){$outdir="./";}
my %config;
###调用程序路径 系统迁移后需相应修改
my $dnadiff="dnadiff";
my $prokka="prokka";
my $blat="blat";
my $fasttree="FastTree";

my @input;
if(defined $dir){
	my @input1=glob "$dir/*.fa";
	my @input2=glob "$dir/*.fasta";
	my @input3=glob "$dir/*.seq";
	@input=(@input1,@input2,@input3);
}else{
	@input=@ARGV;
}
parse_config("$Bin/config.txt",\%config);
if(!defined $identity){$identity=50;}
if(!defined $coverage){$coverage=50;}
if(!defined $type){$type=1;}
if(defined $sum_n){$NUMREPS=int($sum_n/2);$BURNIN=int($sum_n/2);}
if(!defined $sum_n){$sum_n=30000;$NUMREPS=20000;$BURNIN=10000;}
if(!defined $STEP){$STEP=2;}
if($type ne 1 && $type ne 0){$type=1;}

print "Start Analysis\n";
####core gene
my %site;
if(-e $coregene){
	open(F,$coregene);
	while(1){
		my $l=<F>;
		unless($l){last;}
		chomp $l;
		my @a=split",",$l;  ## Core gene分隔符 ,
		foreach  ($a[2]..$a[3]) {
			$site{$a[0]}{$_}++;
		}
	}
	close F;
}else{
	###基因预测
	
	system "$prokka $Ref_s --outdir $outdir/prokka --force --prefix Ref --quiet\n";
	my @seq=@input;
	system "mkdir $outdir/Blat\n";my $strain_num;my %strain;my %core;my %gene;
	foreach my $seq (@seq) {
		unless($seq=~/[0-9a-zA-Z]/){next;}
		$strain_num++;
		my $tn=basename($seq);
		$strain{$tn}++;
		system "$blat $seq $outdir/prokka/Ref.ffn $outdir/Blat/$tn.blat\n";
		open(B,"$outdir/Blat/$tn.blat");
		while(1){
			my $l=<B>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			#print "$a[0],$a[10],$a[9]\n";
			if($a[0]=~/[A-Za-z]/){next;}
			unless($a[0]=~/[0-9]/){next;}
			#print "$a[0],$a[9],$a[8]\n";
			if($a[0]/$a[10]>0.9){
				if(exists $gene{$tn}{$a[9]}){next;}else{$gene{$tn}{$a[9]}=$a[13];$core{$a[9]}++;}
			}
		}
		close B;

	}
	open(F,"$outdir/prokka/Ref.gff");
	open(O,">$outdir/coregene.csv");
	while(1){
			my $l=<F>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			my @b=split"\;",$a[8];
			my $gname=substr($b[0],3);
			unless($gname=~/[0-9a-zA-Z]/){next;}
			#if($l=~/^AE/){print "$l,$a[8],$b[0],$gname,$core{$gname},$strain_num\n";}
			if($core{$gname}>=$strain_num){print O "$a[0],$gname,$a[3],$a[4],$a[6]\n";
			#print  "$a[0],$gname,$a[3],$a[4],$a[6]\n";
				foreach  ($a[3]..$a[4]) {
					
					$site{$a[0]}{$_}++;
				}
			}
	}
	close F;
	close Ol
}

###SNP
my @seq=@input;
#open OUT,">$outdir/snp.sh";
system "mkdir $outdir/mummer\n";my %snp;my %ref;my %chrom;my %csite;
foreach my $seq (@seq) {
       # print OUT "perl $Bin/nucmer_snp.pl $Ref_s $seq $outdir/mummer\n";
		my $name1=basename $seq;
	my $name0=basename $Ref_s;

		my $name=$name1."__".$name0;
		system "$dnadiff $Ref_s $seq -p $outdir/mummer/$name\n";
		open(SNP,"$outdir/mummer/$name.snps")||die;
		while(1){
			my $l=<SNP>;
			unless($l){last;}
			chomp $l;
			my @a=split"\t",$l;
			unless($a[0]=~/[0-9]/){next;}
			unless($a[1]=~/A|T|G|C/){next;}
			my $site=$a[10]."_".$a[0];
			$chrom{$a[10]}++; 
			$csite{$a[0]}++;
			$ref{$site}=$a[1];
			$snp{$name1}{$site}=$a[2];
			#print "$site\n";
		}
		close SNP;

}
close OUT;

###合并SNP

my @strain=keys %snp;
open(O,">$outdir/all_snp.xls");
print O "Ref\tSite\tRef_base\t";
foreach my $strain (@strain) {
	print O "$strain\t";
}
print O "\n";
my @chrom=sort {$a<=>$b}keys %chrom;
my @csite=sort{$a<=>$b} keys %csite;
foreach my $chrom (@chrom) {
	foreach my $csite (@csite) {
	my $site=$chrom."_".$csite;
	unless(exists $ref{$site}){next;}
	print "$site\n";
	print O "$chrom\t$csite\t$ref{$site}\t";
	foreach my $strain (@strain) {
		if(exists $snp{$strain}{$site}){
			print O "$snp{$strain}{$site}\t";
		}else{
			print O "$ref{$site}\t";
		}
	}
	print O "\n";
	}
}
close O;

###合并原有结果
if(defined $old){
	print "Cat old and new result\n";
	open(OLD,"$old");
	my @old;
	my $name=<OLD>;chomp $name;
	my @name=split"\t",$name;
	my $nn=@name;
        $nn=int($nn-3);
		 foreach  (3..$nn) {
			my $site=$name[$_];
			push @old,$site;
		#	print "$site\n";
		 }
		 my %old;
		while (<OLD>) {
				chomp;

				my @a=split"\t",$_;
				if($a[0] eq "Ref"){next;}
				my $nn=@a;
				$nn=int($nn-3);
				#print "$nn\n";
				my $refsite=$a[0]."_".$a[1];
				$old{$refsite}{"Ref"}=$a[2];	
				 foreach  (3..$nn) {
						my $site=$a[$_];
					   $old{$refsite}{$name[$_]}=$site;
						
				}
		}
	close OLD;
	open(OLD,"$outdir/all_snp.xls");
	
	my $name=<OLD>;chomp $name;
	my @name=split"\t",$name;
	my $nn=@name;
        $nn=int($nn-3);
		 foreach  (3..$nn) {
			my $site=$name[$_];
			push @old,$site;
			#print "$site\n";
		 }
		
		while (<OLD>) {
				chomp;

				my @a=split"\t",$_;
				if($a[0] eq "Ref"){next;}
				my $nn=@a;
				$nn=int($nn-3);
				my $refsite=$a[0]."_".$a[1];
				$old{$refsite}{"Ref"}=$a[2];
				
				 foreach  (3..$nn) {
						my $site=$a[$_];
					   $old{$refsite}{$name[$_]}=$site;
						
				}
		}
	close OLD;
	open (NEW,">$outdir/all_snp_new.xls");
	my @site=sort {$a<=>$b}keys %old;
	print NEW "Chrom\tPos\tRef_Pos\t";
	foreach my $old (@old) {
		print NEW "$old\t";
	}
	print NEW "\n";
	my %seq;my %seq_filter;
	foreach my $site  (@site) {
		my @pos=split"_",$site;
		print NEW "$pos[0]\t$pos[1]\t";
		#print "$pos[0]\t$pos[1]\t$site{$pos[0]}{$pos[1]}\n";
		my $ref_pos=$old{$site}{"Ref"};
		$seq{"Ref"}.=$ref_pos;
		if(exists $site{$pos[0]}{$pos[1]}){$seq_filter{"Ref"}.=$ref_pos;}
		print NEW "$ref_pos\t";
		foreach my $old (@old) {
			if(exists $old{$site}{$old}){
				print NEW "$old{$site}{$old}\t";
				$seq{$old}.=$old{$site}{$old};
				
					if(exists $site{$pos[0]}{$pos[1]}){
						$seq_filter{$old}.=$old{$site}{$old};
					}
				
			}else{
				print NEW "$ref_pos\t";
				$seq{$old}.=$ref_pos;
					
					if(exists $site{$pos[0]}{$pos[1]}){
						$seq_filter{$old}.=$ref_pos;
					}
				
			}
		}
		print NEW "\n";
	}
	close NEW;
	open(SEQ,">$outdir/all_snp.fasta");
	my @seqname=keys %seq;
	foreach my $sn (@seqname) {
		print SEQ ">$sn\n$seq{$sn}\n";
	}
	close SEQ;
	open(SEQ,">$outdir/all_snp_filter.fasta");
	my @seqname=keys %seq_filter;
	foreach my $sn (@seqname) {
		print SEQ ">$sn\n$seq_filter{$sn}\n";
	}
	close SEQ;
}else{
	open(OLD,"$outdir/all_snp.xls");
	my @old;my %seq;my %seq_filter;
	my $name=<OLD>;chomp $name;
	my @name=split"\t",$name;
	my $nn=@name;
        $nn=int($nn-3);
		 foreach  (3..$nn) {
			my $site=$name[$_];
			push @old,$site;
			#print "$site\n";
		 }
		 my %old;
		while (<OLD>) {
				chomp;

				my @a=split"\t",$_;
				if($a[0] eq "Ref"){next;}
				my $nn=@a;
				$nn=int($nn-3);
				#print "$nn\n";
				my $refsite=$a[0]."_".$a[1];
				$old{$refsite}{"Ref"}=$a[2];
				$seq{"Ref"}.=$a[2];
				
					if(exists $site{$a[0]}{$a[1]}){
						$seq_filter{"Ref"}.=$a[2];
					}
				
				 foreach  (3..$nn) {
						my $site=$a[$_];
					   $old{$refsite}{$name[$_]}=$site;
					   $seq{$name[$_]}.=$site;
					   
							if(exists $site{$a[0]}{$a[1]}){
								$seq_filter{$name[$_]}.=$site;
							}
					
						
				}
		}
	close OLD;
	
	open(SEQ,">$outdir/all_snp.fasta");
my @seqname=keys %seq;
foreach my $sn (@seqname) {
	print SEQ ">$sn\n$seq{$sn}\n";
}
close SEQ;
open(SEQ,">$outdir/all_snp_filter.fasta");
my @seqname=keys %seq_filter;
foreach my $sn (@seqname) {
	print SEQ ">$sn\n$seq_filter{$sn}\n";
}
close SEQ;
}

my $core_snp_num;
if(defined $old){open(FF,"$outdir/all_snp_new.xls");}else{open(FF,"$outdir/all_snp.xls");}
open (OF,">$outdir/all_snp_filter.xls");
print OF "Ref\tSite\tRef_base\t";
foreach my $strain (@strain) {
	print OF "$strain\t";
}
print OF "\n";
while(1){
	my $line=<FF>;
	unless($line){last;}
	chomp $line;
	if($line=~/^Ref/){next;}
	my @a=split"\t",$line;
	if( exists $site{$a[0]}{$a[1]}){print OF "$line\n";}
	$core_snp_num++;
}
close FF;
print "Core SNP Num : $core_snp_num\n";
print O "Core SNP Num : $core_snp_num\n";
system "mkdir $outdir/result\n";

system "mv $outdir/all_snp.xls $outdir/result/all_snp.xls\n"; #core snp ����
system "mv $outdir/all_snp.fasta $outdir/result/all_snp.xls.fasta\n";
system "mv $outdir/all_snp_filter.xls $outdir/result/\n"; system "mv $outdir/all_snp_filter.fasta $outdir/result/ \n";
system "mv $outdir/coregene.csv $outdir/result/ \n";
if(defined $old){system "mv $outdir/all_snp_new.xls $outdir/result/\n"};

#####构建树
system "$fasttree -quiet <$outdir/result/all_snp_filter.fasta >$outdir/result/all_snp_filter.nwk\n";

print "Finished\n";


#####
my $Time_End= sub_format_datetime(localtime(time()));
print "Running from [$Time_Start] to [$Time_End]\n";

#====================================================================================================================
#  +------------------+
#  |   subprogram     |
#  +------------------+



sub sub_format_datetime #.....
{
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon, $day, $hour, $min, $sec);
}
#####
sub gc{
	my $seq=shift @_;
	my $gc=$seq=~tr/(G|C|g|c)/(G|C|g|c)/;
	my $l=length($seq);


	return ($gc,$l);
}

##parse the software.config file, and check the existence of each software
####################################################
sub parse_config{
	my $conifg_file = shift;
	my $config_p = shift;
	
	my $error_status = 0;
	open IN,$conifg_file || die "fail open: $conifg_file";
	while (<IN>) {
		next if(/#/);
		if (/(\S+)\s*=\s*(\S+)/) {
			my ($software_name,$software_address) = ($1,$2);
			$config_p->{$software_name} = $software_address;

		}
	}
	close IN;
	die "\nExit due to error of software configuration\n" if($error_status);
}
