#Usage:OTU矩阵/Species矩阵 与耐药基因矩阵 比较 寻找相关
if (!require("psych"))
  install.packages("psych") 
 
if (!require("optparse"))
  install.packages("optparse") 

library(psych)
library(optparse)
if (!require("phyloseq"))
  install.packages("phyloseq") 
  library(phyloseq)

#Usage:计算两个矩阵间的关系，用于绘制共现网络
## clean R environment
rm(list = ls())
setwd('./')
## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
## parsing arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input file :colnames is SampleID[Required]"),
  make_option(c("-m", "--meta"), type="character", help="Input file 2:colnames is SampleID [Required]")
 
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$input)) stop('Please input a  table1)')
if(is.null(opts$meta)) stop('Please input a  table2')

matrixfile <- opts$input
mapfile <- opts$meta
otutab = read.delim(matrixfile, row.names=1)
resinfo = read.table(mapfile,row.names=1,header=TRUE)

ps1=phyloseq(otu_table(as.matrix(otutab), taxa_are_rows=TRUE))
ps2=phyloseq(otu_table(as.matrix(resinfo), taxa_are_rows=TRUE))


ps3<-prune_samples(sample_sums(ps1)>50,ps1)#剔除总丰度小于50%的样本
#ps3<-filter_taxa(ps1,function(x) sum (x>0.01)>(0.1*length(x)),TRUE) #仅保留在10%以上的样本中丰度大于0.01的OTU
ps3<-filter_taxa(ps3,function(x) sum (x>0.01)>10,TRUE) #仅保留在5个以上的样本中丰度大于0.01的OTU
ps4<-filter_taxa(ps2,function(x) sum (x>0)>10,TRUE) #仅保留在5个以上的样本中存在的耐药基因

otufilter<-otu_table(ps3)

resfilter<-otu_table(ps4)

CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame(
    from = rownames(cormat)[row(cormat)[ut]],
    to = colnames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}





#两个网络间的关系


filtered_data <- otufilter[,colnames(otufilter) %in% colnames(resfilter)]
filtered_res <- resfilter[,colnames(resfilter) %in% colnames(filtered_data)]


result<-corr.test(t(filtered_data),t(filtered_res),method="spearman",adjust="fdr")
cor_df <- CorrDF(result$r , result$p)
cor_df
#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
#cor_df <- cor_df[which(cor_df$p < 0.05),] # 保留p-value < 0.05 or 0.001的边
write.csv(cor_df,"Psych_result_otu_gene.csv",quote = FALSE,row.names = FALSE)
