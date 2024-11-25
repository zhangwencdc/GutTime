if (!require("ggtree"))
  BiocManager::install("ggtree")
if (!require("treeio"))
  install.packages("treeio") 
if (!require("dplyr"))
  install.packages("dplyr") 
 library("ggtree")
 library("treeio")
 library("dplyr")


 ## parsing arguments
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", help="Input tree file [Required]"),
  make_option(c("-m", "--meta"), type="character", help="Sample Infor file [Required]"),
  
  make_option(c("-o", "--out_dir"), type="character", default='Tree_figure', help="Output directory [default %default]"),
  
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
outpath <- opts$out_dir
# paramenter checking
if(is.null(opts$input_file)) stop('Please input a newick tree file)')
if(is.null(opts$meta)) stop('Please input a sample informatin table(2nd Colname must be "Group"))')
dir.create(outpath)


#input
matrixfile <- opts$input_file
info<-opts$meta
tree<-read.newick(matrixfile,node.label="label")
meta<-read.delim(info,header=TRUE)
# 从文件名提取物种名  
species_name <- sub(".*?\\.([^\\.]+)\\.?.*", "\\1", matrixfile)  
species_name <- gsub("_", " ", species_name) # 替换下划线为空格  

tree_df <- fortify(tree)  

# 合并分组信息  
tree_df <- tree_df %>%  
  left_join(meta, by = c("label" = "ID"))

p1<-  ggtree(tree_df,layout="roundrect") +   
  geom_tree() +  geom_tippoint(aes(label=label,color = Group),size=3)+   scale_color_manual(values = c("P1" = "#7FC97F", "P2" = "#BEAED4", "P3" = "#FDC086","P4"="#FFFF99","P5"="#386CB0","P6"="#F0027F","P7"="#BF5B17"))+
  theme(legend.position = "none")+ annotate("text", x = max(tree$edge.length) * 0.5, y = 5,   label = species_name, size = 2, vjust = -1) 
ggsave(paste(species_name,".pdf", sep="",p1)