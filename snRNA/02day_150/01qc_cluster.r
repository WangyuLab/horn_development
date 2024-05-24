###QC and Cluster ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)

# load dataset, including "barcodes.tsv.gz","features.tsv.gz" and "matrix.mtx.gz" of each data
data <- Read10X("/home/../scRNA/D150horn/data_force_1w")

# creat seurat objects
horn <- CreateSeuratObject(counts = data, project = "day150_horn", min.cells = 10, min.features= 1000)

# calculate mitochondrial gene contenr
mt.genes <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6")
kp=mt.genes %in% rownames(horn)
C<-GetAssayData(object = horn, slot = "counts")
percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
horn <- AddMetaData(horn, percent.mito, col.name = "percent.mito")

# filter cells
VlnPlot(horn, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol= 3)
horn.flit <- subset(horn, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mito < 5 & nCount_RNA < 20000)
VlnPlot(horn.flit, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol= 3)

# remove the ribosome genes
rb.genes <- rownames(horn.flit)[grep("^RP[SL]",rownames(horn.flit))]
length( coi <- rownames( horn.flit )[ which(rowSums(horn.flit@assays$RNA@counts) > 10) ] )
length( coi <- coi[! coi %in% rb.genes ]  )
horn.flit.ribo <- horn.flit[coi,]
horn.flit.ribo <- horn.flit.ribo[! rownames(horn.flit.ribo) %in% rb.genes,]

# normalizing the data
horn.flit.ribo <- NormalizeData(horn.flit.ribo, normalization.method = "LogNormalize", scale.factor = 10000)

# identifiy the highly variable features
horn.flit.ribo <- FindVariableFeatures(horn.flit.ribo, selection.method = "vst", nfeatures = 2000)
var <- VariableFeatures(object = horn.flit.ribo)

# scale the data
all.genes <- rownames(horn.flit.ribo)
horn.flit.ribo <- ScaleData(horn.flit.ribo,features = all.genes,vars.to.regress = "percent.mito")

# romove contamination from the highly variable features
none.genes <- c("HBA","HBG")
var1 <- setdiff(var,mt.genes)
var2 <- setdiff(var1,none.genes)

# perform PCA and determine the dimensionality of the dataset
horn.flit.ribo <- RunPCA(horn.flit.ribo, features = var2)
ElbowPlot(horn.flit.ribo,ndims=50)
DimHeatmap(horn.flit.ribo, dims = 1:25, cells = 500, balanced = TRUE)

# run UMAP and cluster cells
horn.flit.ribo <- RunUMAP(horn.flit.ribo, dims = 1:20)
horn.flit.ribo <- FindNeighbors(horn.flit.ribo, dims = 1:20)
horn.flit.ribo <- FindClusters(horn.flit.ribo, resolution = seq(from =0.1,to = 1.0,by=0.1))

# doublet detect and filter
sweep.res.list <- paramSweep_v3(horn.flit.ribo, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(horn.flit.ribo)*8*1e-6
homotypic.prop <- modelHomotypic(horn.flit.ribo$RNA_snn_res.0.2)
nExp_poi <- round(DoubletRate*ncol(horn.flit.ribo))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
horn.flit.ribo <- doubletFinder_v3(horn.flit.ribo,PCs=1:20,pN=0.25,pK=pK_bcmvn,nExp=nExp_poi.adj,reuse.pANN=F,sct=F)

head(horn.flit.ribo@meta.data)
cells.sub <- subset(horn.flit.ribo@meta.data, DF.classifications_0.25_0.005_665==("Singlet"))
horn.flit.ribo.doublet <- subset(horn.flit.ribo,cells=row.names(cells.sub))

#remove the contaminated subset
sce  <- horn.flit.ribo.doublet
DimPlot(sce,reduction="umap",label=TRUE,group.by="RNA_snn_res.0.3")
Idents(sce)  <- "RNA_snn_res.0.3"
sce <- subset(sce,idents=c(3),inver=TRUE)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
var <- VariableFeatures(object = sce)
all.genes <- rownames(sce)
sce <- ScaleData(sce,features =all.genes,vars.to.regress = "percent.mito")
none.genes <- c("HBA","HBG")
var1 <- setdiff(var,mt.genes)
var2 <- setdiff(var1,none.genes)
sce <- RunPCA(sce, features = var2)
ElbowPlot(sce,ndims=50)
sce <- RunUMAP(sce, dims = 1:20)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = seq(from =0.1,to = 1.0,by=0.1))
DimPlot(sce,label=TRUE,group.by="RNA_snn_res.0.2")

# save image and rds
saveRDS(horn.flit.ribo.doublet, "obj3.rds")
save.image("01qc_cluster.rd")
savehistory("01qc_cluster.r")