#!/usr/bin/perl  
use strict;  
use warnings;  
use File::Basename;  
use File::Path 'make_path';  

# 设置FASTA文件所在的路径  
my $input_dir = $ARGV[0];  
# 设置PlasmidFinder数据库路径  
my $plasmidfinder_db = '/data/zhangwen/bioapp/PlasmidFinder/plasmidfinder_db/';  

# 打开输入目录  
opendir(my $dh, $input_dir) or die "Cannot open directory: $!";  

# 遍历输入目录下的所有FASTA文件  
while (my $filename = readdir($dh)) {  
    if ($filename =~ /\.fasta$/) {  
        # 创建对应的文件夹  
        my $folder_name = basename($filename, '.fasta');  
        my $folder_path = $folder_name;  print "$folder_name\n";
        make_path($folder_path) unless -d $folder_path;  

        # 构建命令  
        my $input_file = "$input_dir/$filename";  
        my $output_dir = $folder_path;  

        # 执行命令  
        my $command = "plasmidfinder.py -i $input_file -o $output_dir -p $plasmidfinder_db -x";  
        system($command) == 0 or warn "System command failed: $!";  
		system "rm -rf $output_dir/tmp\n";
    }  
}  

closedir($dh);