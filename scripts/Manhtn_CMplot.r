#library("ggsci")
library("CMplot")
#https://github.com/YinLiLin/CMplot
#mycolor <- c(pal_d3("category20",alpha = 0.5)(20),pal_d3("category20",alpha = 0.5)(3))
cm<-read.csv("CMP_species.cat.v3",header=TRUE,row.names=1)
#chr_color <- as.data.frame(table(cm$Chromome))


#CMplot(cm[,c("Bacteria","Chromome","Pos","Signals")],plot.type="m",LOG10=TRUE,ylab="-log10(Signals)",type="p",band=1,file="pdf",dpi=300)
CMplot(cm[,c("Bacteria","Chromome","Pos","Signals")],LOG10=TRUE,ylab="-log10(Signals)",type="p",band=1,file="pdf",dpi=300)

library("dplyr")
filtered_df <- cm %>%  
    filter(cm$Signals <= 1e-200)