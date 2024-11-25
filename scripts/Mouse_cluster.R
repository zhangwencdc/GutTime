# 检查BiocManager包，如没有则安装
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
if (!require("meconetcomp"))
  install.packages("meconetcomp") 
library(meconetcomp)
library(microeco)
library(magrittr)
library(ggplot2)

otu = read.delim("Mouse_metaphlan_merge.species", row.names=1)
otu[] <- lapply(otu, function(x) as.numeric(as.character(x))) 
otu[otu <= 0.01] = NA  #仅考虑丰度>0.01%
otu <- otu[(rowSums(is.na(otu)))< 0.9*ncol(otu), ]#仅考虑10%个以上阳性样本的species

#生成函数
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}
#计算OTU间相关性
new<-scale(t(otu))
occor <- corAndPvalue(new, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
taxonomy = read.table("./taxonomy.txt",row.names=1,header=TRUE)

##画整体图
#cols <- c('#00A6FB', '#0582CA', '#fff0f3', '#006494', '#c9e4ca', '#31572c', '#90a955', '#ecf39e', '#a4133c', '#c9184a', '#ff4d6d', '#990033', '#6699cc', '#cc0000', '#00cc66', '#ffcc00', '#666699', '#ff66cc', '#633333', '#99cc00', '#666633', '#339966', '#cc6666', '#6666ff')
cols<-c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499","#44AA99","#999933")


V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
V(igraph)$color <- col_map[V(igraph)$taxon]
V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("All-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("All-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("All-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()

#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"vertices.csv",quote = FALSE,row.names = FALSE)

#net compare
metadata = read.delim("./Mouse_ID-v4.txt",row.names = 1)
#必选，读入otu或宏基因组结果
otutab = read.delim("./Mouse_metaphlan_merge.species", row.names=1)
otutab[] <- lapply(otutab, function(x) as.numeric(as.character(x))) 
# 删除全为0的行  
zero_rows <- apply(otutab, 1, sum) == 0
otutab <- otutab[!zero_rows, ]  

# 删除含有 NA 的行  
otutab <- na.omit(otutab)  

#读入分类关系
taxonomy = read.table("./taxonomy.txt",row.names=1,header=TRUE)
dataset <- microtable$new(sample_table = metadata, otu_table = otutab, tax_table = taxonomy)
dataset_filter <- clone(dataset)
dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.01)

#网络图
t1 <-clone(dataset_filter)
#t1$sample_table %<>%subset(People == "P1")
t1$tidy_dataset()
t1 <- trans_network$new(dataset = t1, cor_method = "spearman", filter_thres = 0.0001)
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy")
t1$save_network(filepath = "AllPeople_network.gexf") #save for Gephi
t1$cal_network_attr()
t1$res_network_attr
t1$get_node_table(node_roles = TRUE)
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
p9<-t1$plot_network()
p9
