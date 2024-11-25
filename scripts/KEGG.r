#KEGG图

#安装包 仅需第一次运行
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#install.packages("tidyjson")
#BiocManager::install("clusterProfiler")

#加载包 每次需运行
library(tidyjson) 
library(jsonlite)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

#加载输入文件
data <- read.table("KEGG_input.txt", sep = "\t", header = T)

#加载kegg注释文件
kegg <- read.csv(file = "kegg_hierarchy.csv")

ggData <- data %>%
  mutate(level3 = kegg$level3[match(ID, kegg$id)],
         level2 = kegg$level2[match(ID, kegg$id)],
         level1 = kegg$level1[match(ID, kegg$id)]) %>%
  
  # 给level1和level2的term排序
  # 注意：只保留table(kegg$level1)里出现的term
  mutate(level1 = factor(level1, 
                          unique(as.character(level1))),
          level2 = factor(level2, 
                          unique(as.character(level2)))
                          ) %>%
  arrange(level1, level2, level3) %>%
  mutate(ID = factor(ID, rev(unique(ID))),
         level3_x = factor(level3, rev(unique(as.character(level3)))))

	 # level1的文字
ggData_l1 <- data.frame(nrow(ggData) - (cumsum(table(ggData$level1)) - table(ggData$level1)/2))
ggData_l1$start <- ggData_l1$Freq - table(ggData$level1)/2
ggData_l1$end <- ggData_l1$Freq + table(ggData$level1)/2

# level2的文字
ggData_l2 <- data.frame(nrow(ggData) - (cumsum(table(ggData$level2)) - table(ggData$level2)/2))
ggData_l2$start <- ggData_l2$Freq - table(ggData$level2)/2
ggData_l2$end <- ggData_l2$Freq + table(ggData$level2)/2


# bar 1
ggplot(ggData) +
  geom_col(mapping = aes(ID, Bloom,fill=level1),
           color = "black", 
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, Bloom, label = Bloom),
            hjust=-0.4, size = 2.5) +
  scale_y_log10(limits = c(1, 1700),
                expand = expansion()) +
 # scale_fill_manual(values = c("#269464", "#CC5618","#2466A2")) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "Bloom") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p1

# bar2
ggplot(ggData) +
  geom_col(mapping = aes(ID, NoBloom,fill=level1),
           color = "black",
           width = 0.75, 
           show.legend = F) +
  geom_text(mapping = aes(ID, NoBloom, label = NoBloom),
            hjust=-0.4, size = 2.5) +
  scale_y_log10(limits = c(1, 1700),
                expand = expansion()) +
#   scale_fill_manual(values = c("#269464", "#CC5618","#2466A2")) +
  coord_flip() + 
  theme_classic() +
  labs(x = NULL, y = NULL, title = "NoBloom") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p2

# text ko
ggplot(ggData) +
  geom_text(mapping = aes(ID, 0, label = ID,color=level1), 
            size = 3, show.legend = F, hjust = 0) +
#  scale_color_manual(values = c("#269464", "#CC5618","#2466A2")) +
  scale_y_continuous(expand = expansion(), limits = c(0,1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Pathway \nIdentifiers") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p0

# text level3
ggplot(ggData) +
  geom_text(mapping = aes(ID, 0, label = level3,color=level1), 
            size = 3, show.legend = F, hjust = 0) +
  #scale_color_manual(values =  c("#269464", "#CC5618","#2466A2")) +
  scale_y_continuous(expand = expansion(), limits = c(0,1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Level 3 of KEGG\nfunctional Category") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p3

# text level2
ggplot(ggData) +
  geom_text(mapping = aes(ID, 0, label = level2,color=level1), 
            size = 3, show.legend = F, hjust = 0) +
  #scale_color_manual(values =  c("#269464", "#CC5618","#2466A2")) +
  scale_y_continuous(expand = expansion(), limits = c(0,1)) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Level 2 of KEGG\nfunctional Category") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p4
# text level1
ggplot(ggData_l1) +
  geom_segment(mapping = aes(x = start+0.1, xend = end-0.1, y = -0.1, yend = -0.1), size = 2)+
  geom_text(mapping = aes(Freq, 0, label = Var1,color=Var1), 
            size = 3, show.legend = F, hjust = 0) +
  scale_y_continuous(expand = expansion(), limits = c(-0.1,1)) +
  scale_x_continuous(expand = expansion()) +
#   scale_fill_manual(values = c("#269464", "#CC5618","#2466A2")) +
  coord_flip() + 
  theme_void() +
  labs(x = NULL, y = NULL, title = "Level 1 of KEGG\nfunctional Category") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 8)) -> p5


cowplot::plot_grid(p0,p1, p2,p3,p4,p5, align = "h", nrow = 1, 
                   rel_widths = c(0.1,0.2,0.2,0.6, 0.5, 0.4))
		   ggsave("KEGGhierarchy.pdf", width = 12, height = 20)