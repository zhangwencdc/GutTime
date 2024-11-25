#https://github.com/YuLab-SMU/MicrobiotaProcess
# install necessary libraries
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}

p <- c("optparse","ggh4x","forcats","gghalves","ggside","ggplot2","treeio","phyloseq","MicrobiotaProcess","ggalluvial","ggtree")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
  BiocManager::install(p)
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')
library(ggplot2)
library(MicrobiotaProcess)
library(ggalluvial)
library(ggtree)
library(gghalves)

## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", help="Input file [Required]"),
  make_option(c("-m", "--meta"), type="character", help="Sample Infor file [Required]"),
  make_option(c("-o", "--out_dir"), type="character", default='MetaFigure', help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
outpath <- opts$out_dir
# paramenter checking
if(is.null(opts$input_file)) stop('Please input a multiple-sample table)')
if(is.null(opts$meta)) stop('Please input a sample informatin table(2nd Colname must be "Group"))')
dir.create(outpath)
# load data
matrixfile <- opts$input_file
meta<-opts$meta

#Input
mpse3 <- mp_import_metaphlan(profile=matrixfile, mapfilename=meta)

#Phylum level冲击图
disbar <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
widforpdf <- ncol(disbar)
p1 <-mpse3 %>%     mp_plot_abundance(        .abundance=Abundance,        .group=Group,          taxa.class = Phylum,         topn = 20,force=TRUE)
ggsave(paste(outpath,"/","Figure1_Phylum_top20.pdf", sep=""), p1, width=ceiling(16+widforpdf/8),height=10, limitsize=FALSE)
#Genus Level冲击图
p2 <-mpse3 %>%     mp_plot_abundance(        .abundance=Abundance,        .group=Group,          taxa.class = Genus,         topn = 20,force=TRUE)
ggsave(paste(outpath,"/","Figure2_Genus_top20.pdf", sep=""), p2, width=ceiling(16+widforpdf/8),height=10, limitsize=FALSE)
#Species/OTU Level冲击图
p3 <-mpse3 %>%     mp_plot_abundance(        .abundance=Abundance,        .group=Group,          taxa.class = OTU,         topn = 20,force=TRUE)
ggsave(paste(outpath,"/","Figure3_Species_top20.pdf", sep=""), p3, width=ceiling(16+widforpdf/8),height=10, limitsize=FALSE)

##alpha多样性
mpse3 %<>% 
mp_cal_alpha(.abundance=Abundance,force=TRUE)
#画图（分组）
f1 <- mpse3 %>% 
      mp_plot_alpha(
        .group=Group, 
        .alpha=c(Observe, Shannon, Simpson, Pielou)
      ) 

ggsave(paste(outpath,"/","Figure4_alpha_diversity_group.pdf", sep=""), f1, width = 5, height = 6)
#画图（不分组）
f2 <- mpse3 %>% mp_plot_alpha(.alpha=c(Observe, Shannon, Simpson, Pielou)) 
ggsave(paste(outpath,"/","Figure5_alpha_diversity.pdf", sep=""), f2, width = 5, height = 6)


##beta多样性
mpse3 %<>%     mp_decostand(.abundance=Abundance)
mpse3 %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
#不分组热力图
b1 <- mpse3 %>% mp_plot_dist(.distmethod = bray)
ggsave(paste(outpath,"/","Figure6_bray_distance_betweenSamples.pdf", sep=""), b1, width = 5, height = 6)
#分组
b2 <- mpse3 %>% mp_plot_dist(.distmethod = bray,.group=Group)
ggsave(paste(outpath,"/","Figure7_bray_distance_betweenSamples_Group.pdf", sep=""), b2, width = 5, height = 6)
#组间比较
b3 <- mpse3 %>% mp_plot_dist(.distmethod = bray, .group = Group, group.test=TRUE, textsize=2)
ggsave(paste(outpath,"/","Figure8_bray_distance_betweenGroup.pdf", sep=""), b3, width = 5, height = 6)
#PCOA分析
mpse3 %<>%     mp_cal_pcoa(.abundance=hellinger, distmethod="bray")# The dimensions of ordination analysis will be added the colData slot (default).
mpse3
c1 <- mpse3 %>%
        mp_plot_ord(
          .ord = pcoa, 
          .group = Group, 
          .color = Group, 
          .size = 1.2,
          .alpha = 1,
          ellipse = TRUE,
     #     show.legend = FALSE # don't display the legend of stat_ellipse
        )
ggsave(paste(outpath,"/","Figure9_PCoA.pdf", sep=""), c1, width = 5, height = 6)

## Sample Cluster:Distance Tree
mpse3 %<>%
       mp_cal_clust(
         .abundance = hellinger, 
         distmethod = "bray",
         hclustmethod = "average", # (UPGAE)
         action = "add" # action is used to control which result will be returned
       )

       sample.clust <- mpse3 %>% mp_extract_internal_attr(name='SampleClust')
sample.clust
p <- ggtree(sample.clust) + 
       geom_tippoint(aes(color=Group)) +
       geom_tiplab(as_ylab = TRUE) +
       ggplot2::scale_x_continuous(expand=c(0, 0.01))
ggsave(paste(outpath,"/","Figure10_Sample_cluster.pdf", sep=""), p, width = 5, height = 6)       

#Taxonmy Tree
mpse3 %<>%
    mp_cal_abundance( # for each samples
      .abundance = Abundance,force=TRUE
    ) %>%
    mp_cal_abundance( # for each groups 
      .abundance=Abundance,
      .group=Group,force=TRUE
    )
      

      mpse3 %<>%
    mp_diff_analysis(
       .abundance = Abundance,
       .group = Group,
       first.test.alpha = 0.05,force=TRUE
    )
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- mpse3 %>% 
               mp_extract_tree(type="taxatree")
taxa.tree
treeout<-paste(outpath,"/Taxonomy.nwk",sep="")
write.tree(as.phylo(taxa.tree),file=treeout)


h1 <- ggtree(
        taxa.tree,
        layout="radial",
        size = 0.3
      ) +
      geom_point(
        data = td_filter(!isTip),
        fill="white",
        size=1,
        shape=21
      )
h2 <- h1 +
      geom_hilight(
        data = td_filter(nodeClass == "Phylum"),
        mapping = aes(node = node, fill = label)
      )
h3<-h2+geom_tiplab(size=2, offset=7.2)

ggsave(paste(outpath,"/","Figure11_TaxonomyTree.pdf", sep=""), h3, width = 5, height = 6)    

####以下是未跑通版本
#
#f.box <- mpse3 %>%
#         mp_plot_diff_boxplot(
#           .group = Group,
#         ) %>%
#         set_diff_boxplot_color(
#           values = c("deepskyblue", "orange"),
#           guide = guide_legend(title=NULL)
#         )
#ggsave(paste(outpath,"/","Figure12_Diff.pdf", sep=""), f.box, width = 5, height = 6)    
#f.bar <- mpse3 %>%
#         mp_plot_diff_boxplot(
#           taxa.class = c(Genus, OTU), # select the taxonomy level to display
#           group.abun = TRUE, # display the mean abundance of each group
#           removeUnknown = TRUE, # whether mask the unknown taxa.
#         ) %>%
#         set_diff_boxplot_color(
#           values = c("deepskyblue", "orange"),
#           guide = guide_legend(title=NULL)
#         )
#ggsave(paste(outpath,"/","Figure13_Diff.pdf", sep=""), f.bar, width = 5, height = 6)    
##aplot::plot_list(f.box, f.bar)