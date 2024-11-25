#!/usr/bin/perl  
use strict;  
use warnings;  

# ���ļ� A �� B  
my $fileA = 's__Paraprevotella_clara_metaphlan.txt';  
my $fileB = 'NC_012920.1-4086';  

# ����һ����ϣ���洢�ļ� A ������  
my %dataA;  

# ��ȡ�ļ� A ������  
open my $fhA, '<', $fileA or die "Could not open file '$fileA': $!";  
while (my $line = <$fhA>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  # �����Ʊ���ָ�  
    my $key = shift @fields;          # ��ȡ��һ����Ϊ��  
    $dataA{$key} = \@fields;          # �洢��������  
}  
close $fhA;  

# ����һ���������洢�ϲ��������  
my @merged_data;  

# ��ȡ�ļ� B �����ݣ����ϲ�  
open my $fhB, '<', $fileB or die "Could not open file '$fileB': $!";  
while (my $line = <$fhB>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  # �����Ʊ���ָ�  
    my $key = shift @fields;          # ��ȡ��һ����Ϊ��  
    
    if (exists $dataA{$key}) {  
        # �ϲ����ݣ������ A ���ҵ���Ӧ�ļ���  
        my @merged_fields = ($key, @{$dataA{$key}}, @fields);  
        push @merged_data, join("\t", @merged_fields);  # ���ϲ�����м�������  
    } else {  
        # ���ֻ�� B ���ҵ���Ӧ�ļ�  
      #  push @merged_data, $line;  
    }  
}  
close $fhB;  

# ���ļ� A ��ʣ�������ӵ��ϲ������  
for my $key (keys %dataA) {  
    unless (grep { $key eq (split /\t/)[0] } @merged_data) {  
        my @merged_fields = ($key, @{$dataA{$key}});  
        push @merged_data, join("\t", @merged_fields);  
    }  
}  

# ��ӡ�ϲ���Ľ��  
foreach my $merged_line (@merged_data) {  
    print "$merged_line\n";  
}