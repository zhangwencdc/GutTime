#renv::activate(project = "~/compare_net")
if (!require("BiocManager"))
  install.packages("BiocManager")
# 检查网络图构建包igraph，如没有则安装
if (!require("igraph"))
  install.packages("igraph")
# 加载网络图构建包igraph
library(igraph)
# 手动安装WGCNA的两个依赖
if (!require("impute"))
  BiocManager::install("impute")
if (!require("preprocessCore"))
  BiocManager::install("preprocessCore")
# 用于计算OTU/ASV之间的相关性
if (!require("WGCNA"))
  install.packages("WGCNA") 
library(WGCNA)

# install necessary libraries
p <- c("optparse","ggplot2","ggsignif","microeco","magrittr","igraph","psych")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p)
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')
library(ggplot2)
library(ggsignif)
library(microeco)

# 使用magrittr包中的管道符
library(magrittr)
library(igraph)
library(ggplot2)
library(psych)
### parsing arguments
#option_list <- list(
#  make_option(c("-i", "--input_file"), type="character", help="Input file [Required]"),
#  make_option(c("-m", "--meta"), type="character", help="Sample Infor file [Required]"),
# # make_option(c("-t", "--taxonomy"), type="character", help="Taxonomy file [Required]"),
#  make_option(c("-o", "--out_dir"), type="character", default='MetaFigure', help="Output directory [default %default]"),
#  make_option(c("-r", "--res"), type="character",  help="Resistance matrix")
#)
#opts <- parse_args(OptionParser(option_list=option_list), args=args)
#outpath <- opts$out_dir
## paramenter checking
#if(is.null(opts$input_file)) stop('Please input a multiple-sample table)')
#if(is.null(opts$meta)) stop('Please input a sample informatin table(2nd Colname must be "Group"))')
#dir.create(outpath)
# load data
matrixfile <- "Time_Humann-cat-genefamilies-cpm.csv.filter"
meta<-"metadata-v3.tsv"


metadata = read.delim(meta,row.names = 1)
otutab = read.delim(matrixfile, row.names=1,header=T)
#taxonomy = read.table(taxo, row.names=1,header=T)

dataset <- microtable$new(sample_table = metadata, otu_table = otutab)
dataset_filter <- clone(dataset)
dataset_filter$filter_taxa(rel_abund = 0.00001)






#beta
d1 <- clone(dataset_filter)
d1$cal_betadiv()
d1$save_betadiv(dirpath="beta_diversity")
t1 <- trans_beta$new(dataset = d1, group = "People", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_PCoA_country.pdf", p7, width =11, height = 5)
t1$cal_ordination(ordination = "PCA")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_PCA_country.pdf", p7, width =11, height = 5)
t1$cal_ordination(ordination = "NMDS")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_NMDS_country.pdf", p7, width =11, height = 5)





#Part4：Network
#spearman
t1 <- trans_network$new(dataset = dataset_filter, cor_method = "spearman", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy")
t1$save_network(filepath = "network.gexf") #save for Gephi
t1$cal_network_attr()
t1$res_network_attr
t1$get_node_table(node_roles = TRUE)
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
p9<-t1$plot_network()
p9
ggsave("P10_net_spearman.pdf", p9, width =11, height = 5)

##分类：Module hubs (模块内部具有高度的连接性, Zi > 2.5 & Pi ≤ 0.62); Connectors (节点主要连接不同模块, Zi ≤ 2.5 & Pi > 0.62); Network hubs (分别具有Module hubs和Connectors的特征, Zi > 2.5 & Pi > 0.62); Peripherals (节点有稀少的连接而且基本是在模块内部, Zi ≤ 2.5 & Pi < 0.62)
#p9<-t1$plot_taxa_roles(use_type = 1)
#ggsave("P9_net.pdf", p9, width =11, height = 5)
#p10<-t1$plot_taxa_roles(use_type = 2)
#ggsave("P10_net.pdf", p10, width =11, height = 5)


###Part5:样本信息 耐药基因与细菌之间的关系
#t2<-trans_env$new(dataset=dataset_filter,add_data=metadata,character2numeric = TRUE) #character2numeric = TRUE将样本信息表中的字节转换成数字，如本身是数字则不应加该参数
#t2$cal_cor(add_abund_table = t1$res_eigen)
#p11<-t2$plot_cor()
#ggsave("P11_heatmap_sampleinfor.pdf", p11, width =11, height = 5)
#
#resinfo = read.table(res,row.names=1,header=TRUE)
#t3<-trans_env$new(dataset=dataset_filter,add_data=t(resinfo))
#t3$cal_cor(add_abund_table = t1$res_eigen)
#p12<-t3$plot_cor()
#ggsave("P11_heatmap_ARG.pdf", p12, width =30, height = 5)
#
##耐药基因和细菌之间的网络关系
##用psych计算相关性 输入的列表行名为样本名
#m1<-merge(otu,res,by="row.names",all=T)
#row.names(m1)<-m1[,1]
#m2<-m1[,-1]
#result<-corr.test(m2,method="spearman",adjust="fdr")
#
#CorrDF <- function(cormat, pmat) {
#  ut <- upper.tri(cormat) 
#  data.frame(
#    from = rownames(cormat)[col(cormat)[ut]],
#    to = rownames(cormat)[row(cormat)[ut]],
#    cor = (cormat)[ut],
#    p = pmat[ut]
#  )
#}
#cor_df <- CorrDF(result$r , result$p)
#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
#cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边
#
