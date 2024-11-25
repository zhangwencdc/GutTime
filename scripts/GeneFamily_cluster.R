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

otu = read.delim("Time_Humann-cat-genefamilies-cpm.csv.filter", row.names=1)
otu[otu <= 100] = NA  #仅考虑丰度大于百万分之100的gene
otu <- otu[(rowSums(is.na(otu)))< 0.5*ncol(otu), ]#仅考虑50%个以上阳性样本的gene

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
#taxonomy = read.table("./taxonomy.txt",row.names=1,header=TRUE)

##画整体图
#cols <- c('#00A6FB', '#0582CA', '#fff0f3', '#006494', '#c9e4ca', '#31572c', '#90a955', '#ecf39e', '#a4133c', '#c9184a', '#ff4d6d', '#990033', '#6699cc', '#cc0000', '#00cc66', '#ffcc00', '#666699', '#ff66cc', '#633333', '#99cc00', '#666633', '#339966', '#cc6666', '#6666ff')
cols<-c("#88CCEE","#CC6677","#DDCC77","#117733","#332288","#AA4499","#44AA99","#999933")


#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
#col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("All-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.frame.color=NA) # 画图
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
#按People分组
metadata = read.delim("./metadata.tsv",row.names = 1)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P1", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
#col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P1-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P1-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P1-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P1_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P1_vertices.csv",quote = FALSE,row.names = FALSE)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P2", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性

cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P2-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P2-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P2-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P2_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P2_vertices.csv",quote = FALSE,row.names = FALSE)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P3", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P3-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P3-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P3-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P3_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P3_vertices.csv",quote = FALSE,row.names = FALSE)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P4", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P4-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P4-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P4-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P4_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P4_vertices.csv",quote = FALSE,row.names = FALSE)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P5", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P5-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P5-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P5-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P5_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P5_vertices.csv",quote = FALSE,row.names = FALSE)


# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P6", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P6-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P6-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P6-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P6_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P6_vertices.csv",quote = FALSE,row.names = FALSE)

# 找到要保留的otu列名
selected_otu_names <- metadata[metadata[,1] == "P7", 0]

# 仅保留otu中列名在selected_otu_names中的行
filtered_otu <- otu[, rownames(selected_otu_names)]
newp<-scale(t(filtered_otu))
occor <- corAndPvalue(newp, use='pairwise', method='spearman') # 计算OTU/ASV之间的spearman相关性
cor_df <- CorrDF(occor$cor , occor$p) # 整理ASV之间的连接关系
cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6),] # 保留spearman相关性绝对值>0.6的边
cor_df <- cor_df[which(cor_df$p < 0.001),] # 保留p-value < 0.001的边

#cor_df <- cor_df[which(abs(cor_df$cor) >= 0.6 & (grepl("g__Prevotella", cor_df$from) | grepl("g__Prevotella", cor_df$to))), ] #仅保留符合特殊字符，比如Prevotella的结果
igraph <- graph_from_data_frame(cor_df, directed= F, vertices= NULL)
length(V(igraph)) # 查看节点数#
length(E(igraph)) # 查看边数
#V(igraph)$taxon <- taxonomy$Phylum[match(V(igraph)$name, rownames(taxonomy))]
col_map <- setNames(cols, unique(V(igraph)$taxon))
# 根据taxonomy表中的Phylum列为节点赋予颜色
#V(igraph)$color <- col_map[V(igraph)$taxon]
#V(igraph)$species <- taxonomy$Species[match(V(igraph)$name, rownames(taxonomy))]
E(igraph)$color[E(igraph)$cor >= 0.6] <- "darkgray" # 正相关则边为深灰色
E(igraph)$color[E(igraph)$cor <= -0.6] <- "red" # 负相关则边为红色
E(igraph)$width <- abs(E(igraph)$cor)*1.2 # 边的粗细与相关系数成正比，进行0.5倍放缩
layout1 <- layout_in_circle(igraph) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(igraph) # fr布局。
layout3 <- layout_on_grid(igraph) # grid布局。
pdf("P7-layout1.pdf", height = 10, width = 10)
plot(igraph, layout=layout1, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P7-layout2.pdf", height = 10, width = 10)
plot(igraph, layout=layout2, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
pdf("P7-layout3.pdf", height = 10, width = 10)
plot(igraph, layout=layout3, vertex.label = V(igraph)$species, vertex.frame.color=NA) # 画图
dev.off()
#保存结果
net.data  <- igraph::as_data_frame(igraph, what = "both")$edges # 提取链接属性
write.csv(net.data,"P7_net.data.csv",quote = FALSE,row.names = FALSE) # 保存结果到本地
vertices  <- igraph::as_data_frame(igraph, what = "both")$vertices # 提取节点属性
write.csv(vertices,"P7_vertices.csv",quote = FALSE,row.names = FALSE)


#net compare
metadata = read.delim("./metadata.tsv",row.names = 1)
#必选，读入otu或宏基因组结果
otutab = read.delim("Time_Humann-cat-genefamilies-cpm.csv.filter", row.names=1)
#读入分类关系
#taxonomy = read.table("./taxonomy.txt",row.names=1,header=TRUE)
dataset <- microtable$new(sample_table = metadata, otu_table = otutab)
dataset_filter <- clone(dataset)
dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.5)

#网络图
t1 <-clone(dataset_filter)

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
ggsave("P10_net_spearman.pdf", p9, width =11, height = 5)
#分People画网络图
###P1
t1 <-clone(dataset_filter)
t1$sample_table %<>%subset(People == "P1")
t1$tidy_dataset()
t1 <- trans_network$new(dataset = t1, cor_method = "spearman", filter_thres = 0.0001)
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy")
t1$save_network(filepath = "P1_network.gexf") #save for Gephi
t1$cal_network_attr()
t1$res_network_attr
t1$get_node_table(node_roles = TRUE)
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
#p9<-t1$plot_network()
#p9
#ggsave("P10_net_spearman_P1.pdf", p9, width =11, height = 5)
###P2
t2 <-clone(dataset_filter)
t2$sample_table %<>%subset(People == "P2")
t2$tidy_dataset()
t2 <- trans_network$new(dataset = t2, cor_method = "spearman", filter_thres = 0.0001)
t2$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t2$cal_module(method = "cluster_fast_greedy")
t2$save_network(filepath = "P2_network.gexf") #save for Gephi
t2$cal_network_attr()
t2$res_network_attr
t2$get_node_table(node_roles = TRUE)
t2$get_edge_table()
# return t2$res_edge_table 
t2$get_adjacency_matrix()
#p9<-t2$plot_network()
#p9
#ggsave("P10_net_spearman_P2.pdf", p9, width =11, height = 5)
###P3
t3 <-clone(dataset_filter)
t3$sample_table %<>%subset(People == "P3")
t3$tidy_dataset()
t3 <- trans_network$new(dataset = t3, cor_method = "spearman", filter_thres = 0.0001)
t3$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t3$cal_module(method = "cluster_fast_greedy")
t3$save_network(filepath = "P3_network.gexf") #save for Gephi
t3$cal_network_attr()
t3$res_network_attr
#t3$get_node_table(node_roles = TRUE)
t3$get_edge_table()
# return t3$res_edge_table 
t3$get_adjacency_matrix()
#p9<-t3$plot_network()
#p9
#ggsave("P10_net_spearman_P3.pdf", p9, width =11, height = 5)
###P4
t4 <-clone(dataset_filter)
t4$sample_table %<>%subset(People == "P4")
t4$tidy_dataset()
t4 <- trans_network$new(dataset = t4, cor_method = "spearman", filter_thres = 0.0001)
t4$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t4$cal_module(method = "cluster_fast_greedy")
t4$save_network(filepath = "P4_network.gexf") #save for Gephi
t4$cal_network_attr()
t4$res_network_attr
t4$get_node_table(node_roles = TRUE)
t4$get_edge_table()
# return t4$res_edge_table 
t4$get_adjacency_matrix()
#p9<-t4$plot_network()
#p9
#ggsave("P10_net_spearman_P4.pdf", p9, width =11, height = 5)
###P5
t5 <-clone(dataset_filter)
t5$sample_table %<>%subset(People == "P5")
t5$tidy_dataset()
t5 <- trans_network$new(dataset = t5, cor_method = "spearman", filter_thres = 0.0001)
t5$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t5$cal_module(method = "cluster_fast_greedy")
t5$save_network(filepath = "P5_network.gexf") #save for Gephi
t5$cal_network_attr()
t5$res_network_attr
t5$get_node_table(node_roles = TRUE)
t5$get_edge_table()
# return t5$res_edge_table 
t5$get_adjacency_matrix()
#p9<-t5$plot_network()
#p9
#ggsave("P10_net_spearman_P5.pdf", p9, width =11, height = 5)
###P6
t6 <-clone(dataset_filter)
t6$sample_table %<>%subset(People == "P6")
t6$tidy_dataset()
t6 <- trans_network$new(dataset = t6, cor_method = "spearman", filter_thres = 0.0001)
t6$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t6$cal_module(method = "cluster_fast_greedy")
t6$save_network(filepath = "P6_network.gexf") #save for Gephi
t6$cal_network_attr()
t6$res_network_attr
t6$get_node_table(node_roles = TRUE)
t6$get_edge_table()
# return t6$res_edge_table 
t6$get_adjacency_matrix()
#p9<-t6$plot_network()
#p9
#ggsave("P10_net_spearman_P6.pdf", p9, width =11, height = 5)
##P7
t7 <-clone(dataset_filter)
t7$sample_table %<>%subset(People == "P7")
t7$tidy_dataset()
t7 <- trans_network$new(dataset = t7, cor_method = "spearman", filter_thres = 0.0001)
t7$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t7$cal_module(method = "cluster_fast_greedy")
t7$save_network(filepath = "P7_network.gexf") #save for Gephi
t7$cal_network_attr()
t7$res_network_attr
t7$get_node_table(node_roles = TRUE)
t7$get_edge_table()
# return t7$res_edge_table 
t7$get_adjacency_matrix()
#p9<-t7$plot_network()
#p9
#ggsave("P10_net_spearman_P7.pdf", p9, width =11, height = 5)



###网络间比较
network<-list()
network$t1<-t1
network$t2<-t2
network$t3<-t3
network$t4<-t4
network$t5<-t5
network$t6<-t6
network$t7<-t7
tmp <- cal_network_attr(network)
write.csv(tmp, "Compare_network_attrubutes.csv")
#比较节点node
tmp <- node_comp(network, property = "name")
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
ggsave("P11_node_overlap_Venn.pdf", g1, width = 7, height = 6)
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard
tmp <- node_comp(network)
tmp1 <- trans_venn$new(tmp)
tmp1$data_summary %<>% .[.[, 1] != 0, ]
g1 <- tmp1$plot_bar()
ggsave("P11_node_overlap_bar.pdf", g1, width = 10, height = 6)
#比较边edge
tmp <- edge_comp(network)
tmp1 <- trans_venn$new(tmp, ratio = "numratio")
g1 <- tmp1$plot_venn(fill_color = FALSE)
ggsave("P11_edge_overlap_Venn.pdf", g1, width = 7, height = 6)
tmp$cal_betadiv(method = "jaccard")
tmp$beta_diversity$jaccard
tmp <- edge_comp(network)
tmp1 <- trans_venn$new(tmp)
tmp1$data_summary %<>% .[.[, 1] != 0, ]#剔除仅在一个网络中存在的边
tmp1$data_summary %<>% .[grepl("&", rownames(.)), ]
g1 <- tmp1$plot_bar()
ggsave("P11_edge_overlap_bar.pdf", g1, width = 10, height = 6)

#Robustness of network
tmp <- robustness$new(network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
    remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)
View(tmp$res_table)
g1<-tmp$plot(linewidth = 1)
ggsave("P11_net_Robustness.pdf", g1, width = 10, height = 6)

#Vulnerability of nodes易损性（vulnerability）：每个节点的易损性为：衡量节点对全局效率（global efficiency）的相对贡献。网络的易损性用网络中节点的最大易损性表示
#Vulnerability：网络中最大的节点脆弱性数值作为该网络的网络脆弱性指标
vul_table <- vulnerability(network)
View(vul_table)
write.csv(vul_table, "Vulnerability_nodes.csv")

#Cohesion
t1 <- cohesionclass$new(network)
View(t1$res_list$sample)
View(t1$res_list$feature)
t1$cal_diff(method = "anova")
g1<-t1$plot(measure = "r_pos")
ggsave("P11_node_Cohesion.pdf", g1, width = 10, height = 6)