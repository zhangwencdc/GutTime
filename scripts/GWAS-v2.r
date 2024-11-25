#conda activate Pysch
if (!require("rMVP"))
  install.packages("rMVP")
  library("rMVP")
if (!require("phyloseq"))
  install.packages("phyloseq") 
  library("phyloseq")
  if (!require("optparse"))
  install.packages("optparse") 

library(optparse)

rm(list = ls())
setwd('./')
## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
## parsing arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input vcf file[Required]"),
  make_option(c("-m", "--meta"), type="character", help="Microbiota data/Phenotype data [Required]")
 
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$input)) stop('Please input a  table1)')
if(is.null(opts$meta)) stop('Please input a  table2')

matrixfile <- opts$input
mapfile <- opts$meta
#matrixfile是vcf文件，mapfile是菌群结构或者样本信息表。注意Sample名必须和vcf文件一致，可能需要手工修改。

#mapfile菌群结构列表进行精简，如为样本信息，可跳过该步
info<-read.table(mapfile,row.names=1,header=TRUE)
info<-t(info)
ps1=phyloseq(otu_table(as.matrix(info), taxa_are_rows=TRUE))
ps2<-filter_taxa(ps1,function(x) sum (x>0.01)>50,TRUE)

n<-nrow(otu_table(ps2))

for (i in 1:n) {
	i
	otu<-t(as.data.frame(otu_table(ps2)[i,]))
	write.csv(otu,file="phe.txt",quote=FALSE)


	#读入mvp
	MVP.Data(fileVCF=matrixfile,filePhe="phe.txt",fileKin=FALSE,filePC=FALSE,out="mvp.vcf",sep.phe=",")
	genotype<-attach.big.matrix("mvp.vcf.geno.desc")
	phenotype<-read.table("mvp.vcf.phe",head=TRUE)
	map<-read.table("mvp.vcf.geno.map",head=TRUE)

	imMVP<-MVP(phe=phenotype,geno=genotype,map=map,threshold=0.05,method=c("GLM","MLM","FarmCPU"))#Significant level: 1.67e-07


}