###  a subcluster analysis of fibroblasts  ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(clustree)

# load the merge data
all1 <- readRDS("/home/../merge/obj5.rds")

# subcluster
fib <- subset(all1,idents=c("FIB"))
fib <- NormalizeData(fib, normalization.method = "LogNormalize", scale.factor = 10000)
fib <- FindVariableFeatures(fib, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fib)
fib <- ScaleData(fib, features = all.genes)
fib <- RunPCA(fib, features = VariableFeatures(object = fib))
ElbowPlot(fib,ndims=50)
fib <- RunHarmony(fib, group.by.vars = "orig.ident", reduction = "pca",dims.use = 1:20, assay.use = "RNA",plot_convergence = TRUE)
fib <- RunUMAP(fib, dims = 1:20, reduction = "harmony", reduction.name = "umap_harmony")
fib <- FindNeighbors(fib,reduction = "harmony", dims = 1:20)
fib <- FindClusters(fib, resolution = seq(from =0.1,to = 1.0,by=0.1))

# annotation
clustree(fib)
Idents(fib) <- "RNA_snn_res.0.3"
levels(fib)
new.cluster.ids=c('FIB1','FIB2','FIB3','FIB4','FIB5','FIB6','FIIB7','FIB8')
names(new.cluster.ids) <- levels(fib)
fib <- RenameIdents(fib, new.cluster.ids)
fib$fib_celltype <- Idents(fib)

# save image and rds
saveRDS(fib, "fib.rds")
save.image("02fib_subcluster.rd")
savehistory("02fib_subcluster.r")



