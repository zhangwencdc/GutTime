#!/usr/bin/perl  
use strict;  
use warnings;  

# 打开文件 A 和 B  
my $fileA = 's__Paraprevotella_clara_metaphlan.txt';  
my $fileB = 'NC_012920.1-4086';  

# 创建一个哈希来存储文件 A 的数据  
my %dataA;  

# 读取文件 A 的内容  
open my $fhA, '<', $fileA or die "Could not open file '$fileA': $!";  
while (my $line = <$fhA>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  # 根据制表符分割  
    my $key = shift @fields;          # 获取第一列作为键  
    $dataA{$key} = \@fields;          # 存储其他数据  
}  
close $fhA;  

# 创建一个数组来存储合并后的数据  
my @merged_data;  

# 读取文件 B 的内容，并合并  
open my $fhB, '<', $fileB or die "Could not open file '$fileB': $!";  
while (my $line = <$fhB>) {  
    chomp $line;  
    my @fields = split /\t/, $line;  # 根据制表符分割  
    my $key = shift @fields;          # 获取第一列作为键  
    
    if (exists $dataA{$key}) {  
        # 合并数据（如果在 A 中找到对应的键）  
        my @merged_fields = ($key, @{$dataA{$key}}, @fields);  
        push @merged_data, join("\t", @merged_fields);  # 将合并后的行加入数组  
    } else {  
        # 如果只在 B 中找到对应的键  
      #  push @merged_data, $line;  
    }  
}  
close $fhB;  

# 将文件 A 中剩余的行添加到合并结果中  
for my $key (keys %dataA) {  
    unless (grep { $key eq (split /\t/)[0] } @merged_data) {  
        my @merged_fields = ($key, @{$dataA{$key}});  
        push @merged_data, join("\t", @merged_fields);  
    }  
}  

# 打印合并后的结果  
foreach my $merged_line (@merged_data) {  
    print "$merged_line\n";  
}