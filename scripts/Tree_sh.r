# 加载必要的库  
library(ggtree)  
library(ggplot2)  
library(cowplot)  
library(ape)  
library(grid)  
library(dplyr)  
library(treeio)  

# 指定包含nwk文件的文件夹  
folder_path <- "C:\\Users\\zhang\\Documents\\Work\\2024Time\\OtherSpecies\\Genome"  
file_list <- list.files(folder_path, pattern = "\\.nwk$", full.names = TRUE)  
meta <- read.delim("C:\\Users\\zhang\\Documents\\Work\\2024Time\\OtherSpecies\\Genome\\ID.txt", header = TRUE)  

# 初始化一个列表以存储所有绘制的图  
plots <- list()  

# 循环遍历所有文件  
for (file in file_list) {  
  # 读取树文件  
  tree <- read.tree(file)  
  # 从文件名提取物种名  
  # 提取文件名  
	file_name <- basename(file)  

	# 使用sub提取物种名  
	species_name <- sub("_(ANI)\\.nwk$", "", file_name)  


  species_name <- gsub("_", " ", species_name) # 替换下划线为空格  

  # 创建树数据框  
  tree_df <- fortify(tree)  

  # 合并分组信息  
  tree_df$label <- sub("\\..*$", "", tree_df$label) #修改label名
 tree_df <- tree_df %>%  
	  left_join(meta, by = c("label" = "ID"))

  # 绘制树  
  p1 <- ggtree(tree_df, layout = "circular") +   
    geom_tree() +  
    geom_tippoint(aes(label = label, color = Group), size = 0.3) +  
    scale_color_manual(values = c("P1" = "#7FC97F", "P2" = "#BEAED4", "P3" = "#FDC086",  
                                   "P4" = "#FFFF99", "P5" = "#386CB0", "P6" = "#F0027F", "P7" = "#BF5B17")) +  
    theme(legend.position = "none") +   
    annotate("text", x = max(tree$edge.length) * 0.5, y = 0,   
             label = species_name, size = 2, vjust = -1)  

  # 保存每个图  
  #ggsave(paste0(species_name, ".pdf"), plot = p1, width =11, height = 5)  

  # 将图添加到列表中  
  plots[[length(plots) + 1]] <- p1  
}  

# 使用cowplot将图合并，每行显示5个图  
combined_plot <- plot_grid(plotlist = plots, ncol = 5)  

# 保存合并后的图形  
ggsave("combined_tree_plot.pdf", plot = combined_plot)