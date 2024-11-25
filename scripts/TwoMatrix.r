#conda activate Pysch
#Usage:OTU����/Species���� ����ҩ������� �Ƚ� Ѱ�����
if (!require("psych"))
  install.packages("psych") 
 
if (!require("optparse"))
  install.packages("optparse") 

library(psych)
library(optparse)
if (!require("phyloseq"))
  install.packages("phyloseq") 
  library(phyloseq)

#Usage:�������������Ĺ�ϵ�����ڻ��ƹ�������
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
if(is.null(opts$meta)) stop('Please input a  table2')#������������ͬ����

matrixfile <- opts$input
mapfile <- opts$meta
otutab = read.delim(matrixfile, row.names=1)

resinfo = read.table(mapfile,row.names=1,header=TRUE)

ps1=phyloseq(otu_table(as.matrix(t(otutab)), taxa_are_rows=TRUE))
ps2=phyloseq(otu_table(as.matrix(t(resinfo)), taxa_are_rows=TRUE))


ps3<-ps1

ps4<-filter_taxa(ps2,function(x) sum (x>0)>3,TRUE) #��������3�����ϵ������д��ڵĻ���

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





#���������Ĺ�ϵ


filtered_data <- otufilter[,colnames(otufilter) %in% colnames(resfilter)]
filtered_res <- resfilter[,colnames(resfilter) %in% colnames(filtered_data)]


result<-corr.test(t(filtered_data),t(filtered_res),method="spearman",adjust="fdr")
cor_df <- CorrDF(result$r , result$p)
cor_df
#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # ����spearman����Ծ���ֵ>0.6�ı�
#cor_df <- cor_df[which(cor_df$p < 0.05),] # ����p-value < 0.05 or 0.001�ı�
write.csv(cor_df,"Psych_result_two_matrix.csv",quote = FALSE,row.names = FALSE)