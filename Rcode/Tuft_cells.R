library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readxl)

d<-readRDS("GSE202100_integrated_epithelial_cells_GEO.RDS")
p0<-DimPlot(d, reduction = "umap", label = T)#+ guides(colour = guide_legend(override.aes = list(size=6)))
p0
RNA<-GetAssay(d, assay = "RNA")
SCT<-GetAssay(d, assay = "SCT")
int<-GetAssay(d, assay = "integrated")

p0<-DimPlot(RNA, reduction = "umap", label = T)+ guides(colour = guide_legend(override.aes = list(size=6)))
p0
p2<-DimPlot(SCT, reduction = "umap", label = T)+ guides(colour = guide_legend(override.aes = list(size=6)))
p2
p3<-DimPlot(int, reduction = "umap", label = T)+ guides(colour = guide_legend(override.aes = list(size=6)))
p3
g<-read_excel("~/Documents/nose/Kotas_JCI insight_paper/Copy of AAP_gene_list_v2.All_genes.with_mouse_orthologs.xlsx")$`Human Gene_Name`

g1<-g[g%in%rownames(d)]
d@meta.data[["cell_type3"]]<-ifelse(d@meta.data[["cell_type2"]]=="Tuft", "Tuft", "Non-Tuft")
d@meta.data[["cell_type3"]]<-factor(d@meta.data[["cell_type3"]])
d@meta.data[["group"]]<-factor(d@meta.data[["group"]])

d@meta.data[["cell_type"]]<-paste0(d@meta.data[["group"]], "_",d@meta.data[["cell_type3"]])

new.cluster.ids <- c("9", "8", "0", "1", "4","2", "3","11", "7","12","5","6", "ionocyte", "13" , "14","tuft" )
names(new.cluster.ids) <- levels(d)
d <- RenameIdents(d, new.cluster.ids)



d1<-subset(d, idents =c("ionocyte","tuft"))

p<-DoHeatmap(d1, features = g1[1:55], group.by = "cell_type",assay = "RNA")
p<-DoHeatmap(d1, features = g1[56:113], group.by = "cell_type",assay = "RNA")

postscript(colormodel="cmyk")
ggsave(filename ='/home/cailu/Documents/nose/Kotas_JCI insight_paper/Figure_4.tiff',p,  width = 60, height =70, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw")
dev.off()


###DElist

dcontrol<-subset(d, group=="Control")
dpolyp<-subset(d, group=="Polyp")


d1control<-subset(d1, group=="Control")
d1polyp<-subset(d1, group=="Polyp")

##perform DE
DE1ctr <- FindMarkers(d1control, ident.1 = "Tuft", group.by="cell_type3", verbose = FALSE)
DE1ctr <-DE1ctr[rownames(DE1ctr)%in%g1,]

DE1polyp <- FindMarkers(d1polyp, ident.1 = "Tuft", group.by="cell_type3", verbose = FALSE)
DE1polyp <-DE1polyp[rownames(DE1polyp)%in%g1,]

##
##perform DE
DEctr <- FindMarkers(dcontrol, ident.1 = "Tuft", group.by="cell_type3", verbose = FALSE)
DEctr <-DEctr[rownames(DEctr)%in%g1,]

DEpolyp <- FindMarkers(dpolyp, ident.1 = "Tuft", group.by="cell_type3", verbose = FALSE)
DEpolyp <-DEpolyp[rownames(DEpolyp)%in%g1,]


##20240305
dx<-subset(d, idents =c("tuft"))
VlnPlot(dx, features="PTGS1",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))

p<-VlnPlot(dx, features="POU2F3",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_POU2F3.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()

#####

p<-VlnPlot(dx, features="SUCNR1",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_SUCNR1.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()
###

p<-VlnPlot(dx, features="CHAT",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_CHAT.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()
####

p<-VlnPlot(dx, features="PTGS1",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_PTGS1.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()

hgpds

####

p<-VlnPlot(dx, features="HPGDS",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_HPGDS.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()

####

p<-VlnPlot(dx, features="BMX",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_BMX.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()


####

p<-VlnPlot(dx, features="GNG13",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_GNG13.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()


####

p<-VlnPlot(dx, features="TRPM5",split.by = "group",split.plot=TRUE)+ xlab("Control       Polyp\nTuft cells")+
  theme(axis.text.x = element_blank(), axis.title = element_text(size=16, face="bold"), axis.text.y=element_text(size=16))


postscript(colormodel="cmyk")
ggsave(filename ='VlnPlot_TRPM5.tiff',p, width =10, height =8, units ="cm",dpi = 600, device='tiff', limitsize = TRUE, compression = "lzw", bg="white")
dev.off()
