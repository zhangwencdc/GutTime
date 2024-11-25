# ���ر�Ҫ�Ŀ�  
library(ggtree)  
library(ggplot2)  
library(cowplot)  
library(ape)  
library(grid)  
library(dplyr)  
library(treeio)  

# ָ������nwk�ļ����ļ���  
folder_path <- "C:\\Users\\zhang\\Documents\\Work\\2024Time\\OtherSpecies\\Genome"  
file_list <- list.files(folder_path, pattern = "\\.nwk$", full.names = TRUE)  
meta <- read.delim("C:\\Users\\zhang\\Documents\\Work\\2024Time\\OtherSpecies\\Genome\\ID.txt", header = TRUE)  

# ��ʼ��һ���б��Դ洢���л��Ƶ�ͼ  
plots <- list()  

# ѭ�����������ļ�  
for (file in file_list) {  
  # ��ȡ���ļ�  
  tree <- read.tree(file)  
  # ���ļ�����ȡ������  
  # ��ȡ�ļ���  
	file_name <- basename(file)  

	# ʹ��sub��ȡ������  
	species_name <- sub("_(ANI)\\.nwk$", "", file_name)  


  species_name <- gsub("_", " ", species_name) # �滻�»���Ϊ�ո�  

  # ���������ݿ�  
  tree_df <- fortify(tree)  

  # �ϲ�������Ϣ  
  tree_df$label <- sub("\\..*$", "", tree_df$label) #�޸�label��
 tree_df <- tree_df %>%  
	  left_join(meta, by = c("label" = "ID"))

  # ������  
  p1 <- ggtree(tree_df, layout = "circular") +   
    geom_tree() +  
    geom_tippoint(aes(label = label, color = Group), size = 0.3) +  
    scale_color_manual(values = c("P1" = "#7FC97F", "P2" = "#BEAED4", "P3" = "#FDC086",  
                                   "P4" = "#FFFF99", "P5" = "#386CB0", "P6" = "#F0027F", "P7" = "#BF5B17")) +  
    theme(legend.position = "none") +   
    annotate("text", x = max(tree$edge.length) * 0.5, y = 0,   
             label = species_name, size = 2, vjust = -1)  

  # ����ÿ��ͼ  
  #ggsave(paste0(species_name, ".pdf"), plot = p1, width =11, height = 5)  

  # ��ͼ��ӵ��б���  
  plots[[length(plots) + 1]] <- p1  
}  

# ʹ��cowplot��ͼ�ϲ���ÿ����ʾ5��ͼ  
combined_plot <- plot_grid(plotlist = plots, ncol = 5)  

# ����ϲ����ͼ��  
ggsave("combined_tree_plot.pdf", plot = combined_plot)