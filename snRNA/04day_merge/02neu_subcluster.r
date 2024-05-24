###  a subcluster analysis of neural cells  ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)

# load the merge data
all1 <- readRDS("/home/../merge/obj5.rds")

# subcluster
neu <- subset(all1,idents=c("NEU"))
neu <- NormalizeData(neu, normalization.method = "LogNormalize", scale.factor = 10000)
neu <- FindVariableFeatures(neu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(neu)
neu <- ScaleData(neu, features = all.genes)

neu <- RunPCA(neu, features = var2)
ElbowPlot(neu,ndims=50)
neu <- RunHarmony(neu, group.by.vars = "orig.ident", reduction = "pca",dims.use = 1:20, assay.use = "RNA",plot_convergence = TRUE)
neu <- RunUMAP(neu, dims = 1:20, reduction = "harmony", reduction.name = "umap_harmony")
neu <- FindNeighbors(neu,reduction = "harmony", dims = 1:20)
neu <- FindClusters(neu, resolution = seq(from =0.1,to = 1.0,by=0.1))

# annotation
clustree(neu)
Idents(neu) <- "RNA_snn_res.0.05"
neu_markers <- FindAllMarkers(neu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
neumarker<- neu_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(neumarker,file="neu_markers.csv")
write.csv(neu_markers,file="neu_allmarkers.csv")

features = c("MKI67","NGFR","SOX10", ##neural crest stem cells marker
	"SNCA","RELN","PAX3",'MPZ', ##neurons and Schwann cells marker
	"POSTN","RUNX2","RXFP2" ## mesenchymal-like cells marker
	)
DotPlot(neu, cols = c('white','#000080'),features = unique(features)) + RotatedAxis()

levels(neu)
new.cluster.ids=c('sub_NEU1','sub_NEU2','Neural_crest')
names(new.cluster.ids) <- levels(neu)
neu <- RenameIdents(neu, new.cluster.ids)
neu$neu_celltype <- Idents(neu)

# save image and rds
saveRDS(neu, "neu.rds")
save.image("02neu_subcluster.rd")
savehistory("02neu_subcluster.r")
