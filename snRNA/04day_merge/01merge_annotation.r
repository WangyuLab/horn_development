### harmony and annotation ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)
library(MySeuratWrappers)

# load dataset
horn1 <- readRDS("/home/../D100/obj2.rds")
horn2 <- readRDS("/home/../D150/obj3.rds")
horn3 <- readRDS("/home/../D200/obj4.rds")

# merge dataset
all <- merge(horn1,y=c(horn2,horn3),add.cell.ids=c("100d","150d","200d"),project="all")

# recalculate the mitochondrial ratio
mt.genes <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6")
kp=mt.genes %in% rownames(all)
C<-GetAssayData(object = all, slot = "counts")
percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
all <- AddMetaData(all, percent.mito, col.name = "percent.mito")
VlnPlot(all, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol= 3)

# normalizing the data, identifiy the highly variable features, scale the data and perform PCA
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
var <- VariableFeatures(object = all)
all.genes <- rownames(all)
all <- ScaleData(all,features =all.genes,vars.to.regress = "percent.mito")
none.genes <- c("HBA","HBG")
var1 <- setdiff(var,mt.genes)
var2 <- setdiff(var1,none.genes)
all <- RunPCA(all, features = var2)
ElbowPlot(all,ndims=50)

# batch correction use Harmony
all1 <- RunHarmony(all, group.by.vars = "orig.ident", reduction = "pca",dims.use = 1:20, assay.use = "RNA",plot_convergence = TRUE)

# dimensionality reduction and cluster cells
all1 <- RunUMAP(all1, dims = 1:20, reduction = "harmony", reduction.name = "umap_harmony")
all1 <- FindNeighbors(all1,reduction = "harmony", dims = 1:20)
all1 <- FindClusters(all1, resolution = seq(from =0.1,to = 1.0,by=0.1))

# choose an appropriate resolution and annotate the cell types
clustree(all1)
DimPlot(all1, reduction = "umap_harmony", group.by = "RNA_snn_res.0.1",label=TRUE)
plot_markers1 <-c("TYR","PROX1","MRC1","RXFP2","PPP2R2B","MYBPC1","MKI67","ITGA7","KDR","KRT5","PDGFRA","OCA2","PAX7","PTPRC")
VlnPlot(all1, features = plot_markers1,stacked=T,pt.size=0)

levels(all1)
new.cluster.ids=c('FIB','FIB','FIB','KRT','END','KRT','MUS1','NEU','MAC','CYC','MEL','LEC','MUS2','KRT','KRT')
names(new.cluster.ids) <- levels(all1)
all1 <- RenameIdents(all1, new.cluster.ids)
all1$merge_celltype <- Idents(all1)

p1 <- DimPlot(all1, reduction = "umap_harmony", group.by = "merge_celltype",label=TRUE)
p2 <- DimPlot(all1, reduction = "umap_harmony", group.by = "merge_celltype")
ggsave(plot=p1,"umap_by_merge_celltype_label.pdf")
ggsave(plot=p2,"umap_by_merge_celltype.pdf")

# plot the proportion of cell types and featureplot of cell markers
Idents(all1) <- "merge_celltype"
tb1 <- table(all1$merge_celltype,all1$orig.ident)
Cellratio1 <- prop.table(tb1,margin = 2)
Cellratio1 <- t(Cellratio1)
Cellratio2 <- prop.table(Cellratio1,margin = 2)
Cellratio2 <- as.data.frame(Cellratio2)
p3 <- ggplot(Cellratio2) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(plot=p3,"mmerge_celltype_proption.pdf")

p4 <- FeaturePlot(all1,features=c("PDGFRA","KRT15","KDR","ITGA7","PPP2R2B","PTPRC","MKI67","TYR","PROX1","MYBPC1"),
	cols = c("lightgrey" ,"#DE1F1F"))
ggsave(plot=p4,"feature_merge_markers.pdf")

# save image and rds
saveRDS(all1, "obj5.rds")
save.image("01merge_annotation.rd")
savehistory("01merge_annotation.r")