###Annotation and Plot ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggplot2)

# load image
load("01qc_cluster.rd")

# choose the appropriate cluster resolution
clustree(all.flit.ribo.doublet)
FeaturePlot(all.flit.ribo.doublet,features=c(
	"PDGFRA", ##fibroblast marker
	"KRT5","KRT15","KRT17", ##keratinocyte marker
	"PECAM1","KDR","CDH5", ##endothelial cell marker
	"PAX7","ITGA7", ## muscle cell marker
	"MKI67", ## prolierating cell marker
	"MYBPC1", ## muscle cell marker
	"PPP2R2B","RELN","SNCA", ## neural cell marker
	"RXFP2", ## horn-specific marker
	"MRC1","PTPRC", ## macrophage marker
	"PROX1","GATA2", ## lymphatic endothelial cell marker
	"OCA2","TYR" ## melanocyte marker
	),col=red)
Idents(all.flit.ribo.doublet) <- "RNA_snn_res.0.6"

# anotate the cell types
levels(all.flit.ribo.doublet)
new.cluster.ids=c('FIB','FIB','FIB','KRT','FIB','FIB','END','MUS1','END','FIB','FIB','CYC',
'MUS2','END','NEU','KRT','RXFP2+_FIB','MAC','CYC','LEC','MEL')
names(new.cluster.ids) <- levels(all.flit.ribo.doublet)
all.ann <- RenameIdents(all.flit.ribo.doublet, new.cluster.ids)
all.ann$celltype <- Idents(all.ann)

# plot 
p1 <- DimPlot(all.ann, reduction = "umap",label=TRUE,group.by="celltype")
p2 <- DimPlot(all.ann, reduction = "umap",group.by="celltype")
p3 <- DimPlot(all.ann, reduction = "umap",label=TRUE,group.by="orig.ident")
p4 <- DimPlot(all.ann, reduction = "umap",group.by="orig.ident")
p5 <- FeaturePlot(all.ann,features=c("PDGFRA","PECAM1","ITGA7","MKI67","MYBPC1","PPP2R2B","RXFP2","MRC1","PROX1","TYR","KDR"),cols = c("lightgrey" ,"#DE1F1F"))
ggsave(plot=p1,"umap_by_celltype_labeled.pdf")
ggsave(plot=p2,"umap_by_celltype.pdf")
ggsave(plot=p3,"umap_by_sample_labeled.pdf")
ggsave(plot=p4,"umap_by_sample.pdf")
ggsave(plot=p5,"featureplot_markers.pdf")

# statistical cell proportion
tb1 <- table(all.ann$celltype,all.ann$orig.ident)
Cellratio <- prop.table(tb1,margin = 2)
Cellratio1 <- t(Cellratio)
Cellratio2 <- prop.table(Cellratio1,margin = 2)
Cellratio2 <- as.data.frame(Cellratio2)
ggplot(Cellratio2) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p6 <- ggplot(Cellratio2) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(plot=p6,"celltype_proption.pdf")

# save image and rds
saveRDS(all.ann, "obj2.rds")
save.image("02annotation_plot.rd")
savehistory("02annotation_plot.r")