The horn bud data of three periods were integrated for analysis

# scATAC data analysis pipeline, mainly with ArchR(v1.0.2)
library('ArchR')
library('Seurat')
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)
addArchRThreads(threads = 16)

# input data from 10x cellranger-atac output
inputFiles <- list.files("{file_pathway}/scATAC-seq/outs",pattern = ".gz",full.names = TRUE)
names(inputFiles) <- gsub("fragments.tsv.gz", "horn", list.files("{file_pathway}/scATAC-seq/outs",pattern = ".gz"))
names(inputFiles)
[1] "100d_horn" "150d_horn" "200d_horn"

#create ArchR object
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = names(inputFiles),filterTSS = 4,filterFrags = 1000,addTileMat = TRUE,addGeneScoreMat = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)

# add doublet score
doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP",LSIMethod = 1)

# create an ArchR object
proj<- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "proj_merge",copyArrows = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)

proj <- filterDoublets(proj)

# QC plot
df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
p1 <- ggPoint(
     x = df[,1],
     y = df[,2],
     colorDensity = TRUE,
     continuousSet = "sambaNight",
     xlabel = "Log10 Unique Fragments",
     ylabel = "TSS Enrichment",
     xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
     ylim = c(0, quantile(df[,2], probs = 0.99))
 ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p2 <- plotFragmentSizes(ArchRProj = proj)
p3 <- plotTSSEnrichment(ArchRProj = proj)
pdf(file="QC_plot.pdf",width = 8, height = 6)
p1+p2+p3
dev.off()

# filter cells
proj <- filterDoublets(ArchRProj = proj)

# dimensional reduction and basic clustering
proj <- addIterativeLSI(ArchRProj = proj,useMatrix = "TileMatrix",name = "IterativeLSI_1",iterations = 2,clusterParams = list(resolution = c(1),sampleCells = 10000,n.start = 10),varFeatures = 25000,dimsToUse = 1:30)
proj <- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI_1",name = "Harmony",groupBy = "Sample")
proj <- addClusters(input = proj,reducedDims = "Harmony",method = "Seurat",name = "Clusters",resolution = 1)

# UMAP embedding
proj <- addUMAP(ArchRProj = proj,reducedDims = "Harmony",name = "UMAP",nNeighbors = 30,minDist = 0.5,metric = "cosine")
p1 <- plotEmbedding(proj,name = "Clusters",embedding = "UMAP",size = 0.7,labelAsFactors=F,labelMeans=F)
p2 <- plotEmbedding(proj,name = "Sample",embedding = "UMAP",size = 0.7,labelAsFactors=F,labelMeans=F)
pdf(file="umap_3sample.pdf",width = 16, height = 6)
p1+p2
dev.off()

#Mapping selected marker genes for different cell types
proj <- addImputeWeights(proj,reducedDims = "Harmony")
markerGenes  <- c("COL1A1","COL3A1","KRT15","KRT17","KRT14","PECAM1","KDR","GRIK2","PPP2R2B","RXFP2","PTPRC","MYH3","MYBPC1"))
p <- plotEmbedding(ArchRProj = proj,colorBy = "GeneScoreMatrix",name = markerGenes, embedding = "UMAP",quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(proj))
pdf(file="horn_3sample_RNAmarker.pdf", width = 8, height = 6)
p
dev.off()

# identification of celltypes based on marker genes
remapClust <- c("C1"="MAC","C2"="FIB","C3"="FIB","C4"="FIB","C5"="FIB","C6"="FIB","C7"="NEU","C8"="FIB","C9"="FIB","C10"="FIB","C11"="FIB","C12"="FIB","C13"="FIB","C14"="FIB","C15"="KRT","C16"="KRT","C17"="FIB","C18"="FIB","C19"="MUS","C20"="MUS","C21"="END","C22"="END")
remapClust <- remapClust[names(remapClust) %in% proj$Clusters]
labelNew <- mapLabels(proj$Clusters, oldLabels = names(remapClust), newLabels = remapClust)
proj$Celltype <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = proj$Clusters)

# plot
p <- plotEmbedding(proj, name = "Celltype", embedding = "UMAP", size = 1.5, labelAsFactors=F, labelMeans=F)
pdf(file="UMAP_3sample_celltype.pdf", width = 8, height = 6)
p1
dev.off()

#Integration of snATAC-seq and snRNA-seq data from horn at three stages.
library(Signac)
library(Seurat)
library(patchwork)
set.seed(1234)

# create combinedATAC object of class Seurat
create_obj <- function(dir) {
count.path <- list.files(path = dir, pattern = "*filtered_peak_bc_matrix.h5", full.names = TRUE)
counts <- Read10X_h5(count.path)
md.path <- list.files(path = dir, pattern = "*singlecell.csv", full.names = TRUE)
md <- read.table(file = md.path, stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
obj <- CreateSeuratObject(counts = counts, assay = "peak", project = 'ATAC',meta.data = md)
return(obj)
}
horn_100 <- create_obj("{file_pathway}/snATAC-seq/outs/horn_100d_outs")
horn_150 <- create_obj("{file_pathway}/snATAC-seq/outs/horn_150d_outs")
horn_200 <- create_obj("{file_pathway}/snATAC-seq/outs/horn_200d_outs")
horn_100$dataset <- 'horn100'
horn_150$dataset <- 'horn150'
horn_200$dataset <- 'horn200'
combinedATAC <- merge(x = horn_100, y = list(horn_150, horn_200), add.cell.ids = c("100", "150", "200"))

# Read fragments files for three periods
allfragments <- CreateFragmentObject(path ="{file_pathway}/scATAC-seq/outs/allfragments.tsv.gz")

# prepare annotation file
library('ensembldb')
gtffile <- "{file_pathway}/Bos_taurus.ARS-UCD1.2.100.gtf"
DB <- ensDbFromGtf(gtf= gtffile)
edb <- EnsDb(DB)
annotations <- genes(edb,filter = ~ gene_biotype == "protein_coding")
genebodyandpromoter.coords <- Extend(x = annotations, upstream = 2000, downstream = 0)
genebodyandpromoter.coords$gene_name[is.na(genebodyandpromoter.coords$gene_name)] <- genebodyandpromoter.coords$gene_id; genebodyandpromoter.coords
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)

# Create a gene activity matrix
gene.activities <- FeatureMatrix(fragments = allfragments,features = genebodyandpromoter.coords,cells = colnames(combinedATAC))
rownames(gene.activities) <- gene.key[rownames(gene.activities)]
combinedATAC[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

#
combinedATAC$tech <- "atac"
DefaultAssay(combinedATAC) <- "ACTIVITY"
combinedATAC <- FindVariableFeatures(combinedATAC)
combinedATAC <- NormalizeData(combinedATAC)
combinedATAC <- ScaleData(combinedATAC)
DefaultAssay(combinedATAC) <- "peak"
VariableFeatures(combinedATAC) <- names(which(Matrix::rowSums(combinedATAC) > 100))
combinedATAC <- RunTFIDF(combinedATAC)
combinedATAC <- FindTopFeatures(combinedATAC, min.cutoff = 'q0')
combinedATAC <- RunSVD(combinedATAC)

# Read the snRNA-seq data previously generated by Seurat
horn.rna <- readRDS("{file_pathway}/snRNA-seq/horn/merge_celltype.rds")
horn.rna$tech <- "rna"

#Identify anchors
transfer.anchors <- FindTransferAnchors(reference = horn.rna,query = combinedATAC,features = VariableFeatures(object = horn.rna),reference.assay = "RNA",query.assay = "ACTIVITY",reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = horn.rna$merge_celltype,weight.reduction = combinedATAC[["lsi"]],dims = 2:50)
combinedATAC <- AddMetaData(combinedATAC, metadata = celltype.predictions)
pdf(file="prediction.score.pdf", width = 8, height = 6)
hist(combinedATAC$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()
horn.atac.filtered <- subset(combinedATAC,subset = prediction.score.max > 0.5)
p1 <- DimPlot(horn.atac.filtered,group.by = "predicted.id",label = TRUE, repel = TRUE)+ggtitle("scATAC-seq cells") + NoLegend()
p2 <- DimPlot(horn.rna, group.by = "merge_celltype", label = TRUE, repel = TRUE) +ggtitle("scRNA-seq cells") +NoLegend()
pdf(file="atac_rna_celltype_match.pdf", width = 12, height = 6)
p1+p2
dev.off() 

#Co-embedding snATAC-seq and snRNA-seq datasets
genes.use <- VariableFeatures(horn.rna)
refdata <- GetAssayData(horn.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors,refdata = refdata,weight.reduction =combinedATAC[["lsi"]],dims=2:50)
combinedATAC[["RNA"]] <- imputation
coembed <- merge(x = horn.rna, y = combinedATAC)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$cellType <- ifelse(!is.na(coembed$cellType), coembed$cellType, coembed$predicted.id)
p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "cellType", label = TRUE, repel = TRUE)
pdf(file="atac_rna_coembed.pdf", width = 12, height = 6)
p1+p2
dev.off()
p3 <- DimPlot(coembed, reduction="umap", split.by = "tech", group.by = "cellType", label = TRUE, repel = TRUE) + NoLegend()
pdf(file="atac_rna_coembed_split.pdf", width = 10, height = 6)
p3
dev.off()

#FIB cells were extracted
idxFB <- BiocGenerics::which(proj$Celltype %in% c("FIB"))
cellsSample <- proj$cellNames[idxFB]
proj[cellsSample, ]
saveArchRProject(ArchRProj = proj[cellsSample, ],outputDirectory = './save_proj_FIBsub',load = FALSE,overwrite = TRUE)
proj <- addIterativeLSI(ArchRProj = proj,useMatrix = "TileMatrix",name = "IterativeLSI_0.5",iterations = 2,clusterParams = list(resolution = c(0.5),sampleCells = 10000,n.start = 10),varFeatures = 25000,dimsToUse = 1:30,force = TRUE)
proj <- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI_0.5",name = "Harmony",groupBy = "Sample",force = TRUE)
proj <- addClusters(input = proj,reducedDims = "Harmony",method = "Seurat",name = "Clusters",resolution = 0.5,force = TRUE)
proj <- addUMAP(ArchRProj = proj,reducedDims = "Harmony",name = "UMAP",nNeighbors = 30,minDist =0.5,metric = "cosine",force = TRUE)
pdf(file="UMAP_FIB.pdf", width = 8, height = 6)
p1
dev.off()

#The corresponding FIB subclusters are extracted
#RXFP2 upstream specifically open peak
gr <- GRanges(seqnames = "chr12", ranges = IRanges(start =c(28774746),end=c(29774747)))
p <- plotBrowserTrack(ArchRProj = proj, groupBy = "cellType",region=gr,geneSymbol = markerGenes,loops = getCoAccessibility(proj))
pdf("CoAcc_RXFP2.pdf")
grid::grid.draw(p)
dev.off()

#TF motif scanning of specific open peak
bedtools getfasta -fi Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa -bed usepeak.bed -fo usepeak.fasta
mast JASPAR2022_CORE_vertebrates_redundant.meme usepeak.fasta -nostatus -hit_list -best > usepeak.result





















































































