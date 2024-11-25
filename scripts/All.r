#install.packages("meconetcomp")
#renv::activate(project = "~/compare_net")
#install.packages("meconetcomp")
#renv::activate(project = "~/compare_net")
#install.packages("meconetcomp")
#install.packages("networkD3")
#install.packages("ggraph")
#install.packages("rgexf")
## heatmap使用
#install.packages("pheatmap")
## microeco包的依赖包之一
#install.packages("aplot")
##用来进行 Duncan's new multiple range test
#install.packages("agricolae")
#install.packages("devtools")
#devtools::install_github
#renv::snapshot()
#以下资料来源：https://chiliubio.github.io/microeco_tutorial 根据自有数据修改调整
#
library(microeco)
library(meconetcomp)
# 使用magrittr包中的管道符
library(magrittr)
library(igraph)
library(ggplot2)

metadata = read.delim("./metadata.tsv",row.names = 1)
#必选，读入otu或宏基因组结果
otutab = read.delim("./otutab.txt", row.names=1)
#读入分类关系
taxonomy = read.table("./taxonomy.txt",row.names=1,header=TRUE)
dataset <- microtable$new(sample_table = metadata, otu_table = otutab, tax_table = taxonomy)
dataset_filter <- clone(dataset)
dataset_filter$filter_taxa(rel_abund = 0.0001, freq = 0.1)

# Part1：Phylum Abundance  Top10 如需要其他级别，直接替换成Genus
t1 <- trans_abund$new(dataset = dataset_filter, taxrank = "Phylum", ntaxa = 10)
p1<-t1$plot_bar(others_color = "grey70", facet = "People", xtext_keep = FALSE, legend_text_italic = FALSE)
p1
ggsave("P1_Phylum_abundance.pdf", p1, width =11, height = 5)

t1 <- trans_abund$new(dataset = dataset_filter, taxrank = "Phylum", ntaxa = 10, groupmean = "People")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
p2<-g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
ggsave("P2_Phylum_abundance_People.pdf", p2, width =11, height = 5)

t1 <- trans_abund$new(dataset = dataset_filter, taxrank = "Phylum", ntaxa = 10)
p3<-t1$plot_box(group = "People", xtext_angle = 30)
p3
ggsave("P3_Phylum_boxplot.pdf", p3, width =11, height = 5)

t1 <- trans_abund$new(dataset = dataset_filter, taxrank = "Phylum", ntaxa = 10)
g1 <- t1$plot_heatmap(facet = "People", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
g1
p4<-g1 + theme(axis.text.y = element_text(face = 'italic'))
p4
ggsave("P4_Phylum_heatmap.pdf", p4, width =11, height = 5)

dataset1 <- dataset_filter$merge_samples(use_group = "People")
t1 <- trans_venn$new(dataset = dataset1)


# Part2：alpha
t1 <- trans_alpha$new(dataset = dataset_filter, group = "People")

t1$cal_diff(method = "wilcox")
p5<-t1$plot_alpha(measure = "Shannon",add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE)
ggsave("P5_Shannon.pdf", p5, width =11, height = 5)
p6<-t1$plot_alpha(measure = "Pielou",add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE)
ggsave("P6_Pielou.pdf", p6, width =11, height = 5)

#beta
d1 <- clone(dataset_filter)
d1$cal_betadiv()
t1 <- trans_beta$new(dataset = d1, group = "People", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_PCoA.pdf", p7, width =11, height = 5)
t1$cal_ordination(ordination = "PCA")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_PCA.pdf", p7, width =11, height = 5)
t1$cal_ordination(ordination = "NMDS")
class(t1$res_ordination)
p7<-t1$plot_ordination(plot_color = "People", plot_shape = "People", plot_type = c("point", "ellipse"))
p7
ggsave("P7_NMDS.pdf", p7, width =11, height = 5)



#Part3：Lefse
t1 <- trans_diff$new(dataset = dataset_filter, method = "lefse", group = "People", alpha = 0.01, lefse_subgroup = NULL)
p7<-t1$plot_diff_bar(threshold = 4)#输出LDA>4的结果
p8<-t1$plot_diff_abund(use_number = 1:30)#输出top30的结果
ggsave("P8_LEFSE.pdf", p7, width =11, height = 5)
ggsave("P8_LEFSE_top30.pdf", p8, width =11, height = 5)
#制定Genus level比较
t1 <- trans_diff$new(dataset = dataset_filter, method = "wilcox", group = "People", taxa_level = "Genus", filter_thres = 0.001)
t1$res_diff %<>% subset(Significance %in% "***")
t1$plot_diff_abund(use_number = 1:10, add_sig = T, add_sig_label = "Significance")

#Part4：Network
#spearman
t1 <- trans_network$new(dataset = dataset_filter, cor_method = "spearman", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
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

#sparcc
#NetComi包安装有问题
#t1 <- trans_network$new(dataset = dataset_filter, cor_method = "sparcc", use_sparcc_method = "NetCoMi", filter_thres = 0.001)
#t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
#t1$cal_module(method = "cluster_fast_greedy")
#t1$save_network(filepath = "network.gexf") #save for Gephi
#t1$cal_network_attr()
#t1$res_network_attr
#t1$get_node_table(node_roles = TRUE)
#t1$get_edge_table()
## return t1$res_edge_table 
#t1$get_adjacency_matrix()
#p9<-t1$plot_network()
#p9
#ggsave("P10_net_spearman.pdf", p9, width =11, height = 5)


#分类：Module hubs (模块内部具有高度的连接性, Zi > 2.5 & Pi ≤ 0.62); Connectors (节点主要连接不同模块, Zi ≤ 2.5 & Pi > 0.62); Network hubs (分别具有Module hubs和Connectors的特征, Zi > 2.5 & Pi > 0.62); Peripherals (节点有稀少的连接而且基本是在模块内部, Zi ≤ 2.5 & Pi < 0.62)
p9<-t1$plot_taxa_roles(use_type = 1)
ggsave("P9_net.pdf", p9, width =11, height = 5)
p10<-t1$plot_taxa_roles(use_type = 2)
ggsave("P10_net.pdf", p10, width =11, height = 5)

#分People画网络图
###P1
t1 <-clone(dataset_filter)
t1$sample_table %<>%subset(People == "P1")
t1$tidy_dataset()
t1 <- trans_network$new(dataset = t1, cor_method = "spearman", filter_thres = 0.001)
t1$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy")
t1$save_network(filepath = "P1_network.gexf") #save for Gephi
t1$cal_network_attr()
t1$res_network_attr
t1$get_node_table(node_roles = TRUE)
t1$get_edge_table()
# return t1$res_edge_table 
t1$get_adjacency_matrix()
p9<-t1$plot_network()
p9
ggsave("P10_net_spearman_P1.pdf", p9, width =11, height = 5)
###P2
t2 <-clone(dataset_filter)
t2$sample_table %<>%subset(People == "P2")
t2$tidy_dataset()
t2 <- trans_network$new(dataset = t2, cor_method = "spearman", filter_thres = 0.001)
t2$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t2$cal_module(method = "cluster_fast_greedy")
t2$save_network(filepath = "P2_network.gexf") #save for Gephi
t2$cal_network_attr()
t2$res_network_attr
t2$get_node_table(node_roles = TRUE)
t2$get_edge_table()
# return t2$res_edge_table 
t2$get_adjacency_matrix()
p9<-t2$plot_network()
p9
ggsave("P10_net_spearman_P2.pdf", p9, width =11, height = 5)
###P3
t3 <-clone(dataset_filter)
t3$sample_table %<>%subset(People == "P3")
t3$tidy_dataset()
t3 <- trans_network$new(dataset = t3, cor_method = "spearman", filter_thres = 0.001)
t3$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t3$cal_module(method = "cluster_fast_greedy")
t3$save_network(filepath = "P3_network.gexf") #save for Gephi
t3$cal_network_attr()
t3$res_network_attr
t3$get_node_table(node_roles = TRUE)
t3$get_edge_table()
# return t3$res_edge_table 
t3$get_adjacency_matrix()
p9<-t3$plot_network()
p9
ggsave("P10_net_spearman_P3.pdf", p9, width =11, height = 5)
###P4
t4 <-clone(dataset_filter)
t4$sample_table %<>%subset(People == "P4")
t4$tidy_dataset()
t4 <- trans_network$new(dataset = t4, cor_method = "spearman", filter_thres = 0.001)
t4$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t4$cal_module(method = "cluster_fast_greedy")
t4$save_network(filepath = "P4_network.gexf") #save for Gephi
t4$cal_network_attr()
t4$res_network_attr
t4$get_node_table(node_roles = TRUE)
t4$get_edge_table()
# return t4$res_edge_table 
t4$get_adjacency_matrix()
p9<-t4$plot_network()
p9
ggsave("P10_net_spearman_P4.pdf", p9, width =11, height = 5)
###P5
t5 <-clone(dataset_filter)
t5$sample_table %<>%subset(People == "P5")
t5$tidy_dataset()
t5 <- trans_network$new(dataset = t5, cor_method = "spearman", filter_thres = 0.001)
t5$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t5$cal_module(method = "cluster_fast_greedy")
t5$save_network(filepath = "P5_network.gexf") #save for Gephi
t5$cal_network_attr()
t5$res_network_attr
t5$get_node_table(node_roles = TRUE)
t5$get_edge_table()
# return t5$res_edge_table 
t5$get_adjacency_matrix()
p9<-t5$plot_network()
p9
ggsave("P10_net_spearman_P5.pdf", p9, width =11, height = 5)
###P6
t6 <-clone(dataset_filter)
t6$sample_table %<>%subset(People == "P6")
t6$tidy_dataset()
t6 <- trans_network$new(dataset = t6, cor_method = "spearman", filter_thres = 0.001)
t6$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t6$cal_module(method = "cluster_fast_greedy")
t6$save_network(filepath = "P6_network.gexf") #save for Gephi
t6$cal_network_attr()
t6$res_network_attr
t6$get_node_table(node_roles = TRUE)
t6$get_edge_table()
# return t6$res_edge_table 
t6$get_adjacency_matrix()
p9<-t6$plot_network()
p9
ggsave("P10_net_spearman_P6.pdf", p9, width =11, height = 5)
##P7
t7 <-clone(dataset_filter)
t7$sample_table %<>%subset(People == "P7")
t7$tidy_dataset()
t7 <- trans_network$new(dataset = t7, cor_method = "spearman", filter_thres = 0.001)
t7$cal_network(COR_p_thres = 0.01, COR_cut = 0.6)
t7$cal_module(method = "cluster_fast_greedy")
t7$save_network(filepath = "P7_network.gexf") #save for Gephi
t7$cal_network_attr()
t7$res_network_attr
t7$get_node_table(node_roles = TRUE)
t7$get_edge_table()
# return t7$res_edge_table 
t7$get_adjacency_matrix()
p9<-t7$plot_network()
p9
ggsave("P10_net_spearman_P7.pdf", p9, width =11, height = 5)

#网络图参数说明
#vertex.label  = NA 不加label
#V(t1$res_network)$species<- taxonomy$Species[match(V(t1$res_network)$name, rownames(taxonomy))]
#t1$plot_network(vertex.label = V(t1$res_network)$species) 改lable名
#t1$plot_network(vertex.label = V(t1$res_network)$species,layout=layout_in_circle) 修改排布方式

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
write.csv(tmp, "network_attrubutes.csv")
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


##Part5:样本信息 耐药基因与细菌之间的关系
t2<-trans_env$new(dataset=dataset_filter,add_data=metadata,character2numeric = TRUE) #character2numeric = TRUE将样本信息表中的字节转换成数字，如本身是数字则不应加该参数
t2$cal_cor(add_abund_table = t1$res_eigen)
p11<-t2$plot_cor()
ggsave("P11_heatmap_sampleinfor.pdf", p11, width =11, height = 5)

resinfo = read.table("ARG.stat.detail.matrix",row.names=1,header=TRUE)
t3<-trans_env$new(dataset=dataset_filter,add_data=t(resinfo))
t3$cal_cor(add_abund_table = t1$res_eigen)
p12<-t3$plot_cor()
ggsave("P11_heatmap_ARG.pdf", p12, width =30, height = 5)

#用psych计算相关性
result<-corr.test(otutab,method="spearman",adjust="fdr")