###QC and Cluster ###

# load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(DoubletFinder)


# load dataset, including "barcodes.tsv.gz","features.tsv.gz" and "matrix.mtx.gz" of each data
horn.data <- Read10X("/home/../D100horn/data/")
skin.data <- Read10X("/home/../D100skin/data/")

# creat seurat objects
horn <- CreateSeuratObject(counts = horn.data, project = "day100_horn", min.cells = 5, min.features= 1000)
skin <- CreateSeuratObject(counts = skin.data, project = "day100_skin", min.cells = 5, min.features= 1000)

# merging data
all <- merge(horn,skin)

# calculate mitochondrial gene contenr
mt.genes <- c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6")
kp=mt.genes %in% rownames(all)
C<-GetAssayData(object = all, slot = "counts")
percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
all <- AddMetaData(all, percent.mito, col.name = "percent.mito")

# filter cells
VlnPlot(all, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol= 3)
all.flit <- subset(all, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mito < 5 & nCount_RNA < 20000)
VlnPlot(all.flit, features = c("nFeature_RNA","nCount_RNA", "percent.mito"), ncol= 3)

# remove the ribosome genes
rb.genes <- rownames(all.flit)[grep("^RP[SL]",rownames(all.flit))]
length( coi <- rownames( all.flit )[ which(rowSums(all.flit@assays$RNA@counts) > 10) ] )
length( coi <- coi[! coi %in% rb.genes ]  )
all.flit.ribo <- all.flit[coi,]
all.flit.ribo <- all.flit.ribo[! rownames(all.flit.ribo) %in% rb.genes,]

# normalizing the data
all.flit.ribo <- NormalizeData(all.flit.ribo, normalization.method = "LogNormalize", scale.factor = 10000)

# identifiy the highly variable features
all.flit.ribo <- FindVariableFeatures(all.flit.ribo, selection.method = "vst", nfeatures = 2000)
var <- VariableFeatures(object = all.flit.ribo)

# scale the data
all.genes <- rownames(all.flit.ribo)
all.flit.ribo <- ScaleData(all.flit.ribo,features =all.genes,vars.to.regress = "percent.mito")

# romove contamination from the highly variable features
none.genes <- c("HBA","HBG")
var1 <- setdiff(var,mt.genes)
var2 <- setdiff(var1,none.genes)

# perform PCA and determine the dimensionality of the dataset
all.flit.ribo <- RunPCA(all.flit.ribo, features = var2)
ElbowPlot(all.flit.ribo,ndims=50)
DimHeatmap(all.flit.ribo, dims = 1:25, cells = 500, balanced = TRUE)

# run UMAP and cluster cells
all.flit.ribo <- RunUMAP(all.flit.ribo, dims = 1:20)
all.flit.ribo <- FindNeighbors(all.flit.ribo, dims = 1:20)
all.flit.ribo <- FindClusters(all.flit.ribo, resolution = seq(from =0.1,to = 1.0,by=0.1))

# doublet detect and filter
sweep.res.list <- paramSweep_v3(all.flit.ribo, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(all.flit.ribo)*8*1e-6
homotypic.prop <- modelHomotypic(all.flit.ribo$RNA_snn_res.0.2)
nExp_poi <- round(DoubletRate*ncol(all.flit.ribo))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
all.flit.ribo <- doubletFinder_v3(all.flit.ribo,PCs=1:20,pN=0.25,pK=pK_bcmvn,nExp=nExp_poi.adj,reuse.pANN=F,sct=F)

head(all.flit.ribo@meta.data)
cells.sub <- subset(all.flit.ribo@meta.data,DF.classifications_0.25_0.005_1835==("Singlet"))
all.flit.ribo.doublet <- subset(all.flit.ribo,cells=row.names(cells.sub))

# save image and rds
saveRDS(all.flit.ribo.doublet, "obj1.rds")
save.image("01qc_cluster.rd")
savehistory("01qc_cluster.r")