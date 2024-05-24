
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
inputFiles <- list.files("{file_pathway}/horn_vs_skin", pattern = ".gz",full.names = TRUE)
names(inputFiles) <- gsub(".fragments.tsv.gz", "", list.files("{file_pathway}/horn_vs_skin",pattern = ".gz"))

# create ArchR object
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = names(inputFiles),filterTSS = 4,filterFrags = 1000,addTileMat = TRUE,addGeneScoreMat = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)

# add doublet score
doubScores <- addDoubletScores(input = ArrowFiles, k = 10,knnMethod = "UMAP",LSIMethod = 1)

# create an ArchR object
proj<- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "proj",copyArrows = TRUE,geneAnnotation = geneAnnotation,genomeAnnotation = genomeAnnotation)

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
proj <- addIterativeLSI(ArchRProj = proj,useMatrix = "TileMatrix",name = "IterativeLSI_0.8",iterations = 2,clusterParams = list(resolution = c(0.8),sampleCells = 10000,n.start = 10),varFeatures = 25000,dimsToUse = 1:30)
proj <- addHarmony(ArchRProj = proj,reducedDims = "IterativeLSI",name = "Harmony",groupBy = "Sample")
proj <- addClusters(input = proj,reducedDims = "IterativeLSI_0.8",method = "Seurat",name = "Clusters",resolution = 0.8)

# UMAP embedding
proj <- addUMAP(ArchRProj = proj,reducedDims = "Harmony",name = "UMAP",nNeighbors = 30,minDist = 0.5,metric = "cosine",force = TRUE)
p1 <- plotEmbedding(proj,name = "Clusters",embedding = "UMAP",size = 0.7,labelAsFactors=F,labelMeans=F)
p2 <- plotEmbedding(proj,name = "Sample",embedding = "UMAP",size = 0.7,labelAsFactors=F,labelMeans=F)
pdf(file="umap.pdf",width = 16, height = 6)
p1+p2
dev.off()
 
# Mapping selected marker genes for different cell types
proj <- addImputeWeights(proj)
markerGenes  <- c("COL1A1","COL3A1","KRT15","KRT17","KRT14","PECAM1","KDR","GRIK2","PPP2R2B","RXFP2","PTPRC","MYH3","MYBPC1")
p <- plotEmbedding(ArchRProj = proj,colorBy = "GeneScoreMatrix",name = markerGenes, embedding = "UMAP",quantCut = c(0.01, 0.95),imputeWeights = getImputeWeights(proj))
pdf(file="horn_vs_skin_genemarker.pdf", width = 8, height = 6)
p
dev.off()

# identification of celltypes based on marker genes
remapClust <- c("C1"="MAC","C2"="FIB","C3"="FIB","C4"="FIB","C5"="FIB","C6"="FIB","C7"="RXFP2+_FIB","C8"="FIB","C9"="FIB","C10"="NEU","C11"="FIB","C12"="FIB","C13"="MUS","C14"="MUS","C15"="MUS","C16"="KRT","C17"="END")
remapClust <- remapClust[names(remapClust) %in% proj$Clusters]
labelNew <- mapLabels(proj$Clusters, oldLabels = names(remapClust), newLabels = remapClust)
proj$Celltype <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = proj$Clusters)

#plot
p <- plotEmbedding(proj, name = "Celltype", embedding = "UMAP", size = 1.5, labelAsFactors=F, labelMeans=F)
pdf(file="UMAP_celltype.pdf", width = 8, height = 6)
p1
dev.off()

# peak calling with macs2 v2.2.7
pathToMacs2 <- findMacs2()
proj_peakclu <- addGroupCoverages(ArchRProj = proj, groupBy = "Celltype")
proj <- addReproduciblePeakSet(ArchRProj = proj,groupBy = "Celltype",pathToMacs2 = pathToMacs2,genomeSize = 2.7e+09)
proj <- addPeakMatrix(proj)

#celltype specific peaks
markersPeaks <- getMarkerFeatures(ArchRProj = proj,useMatrix = "PeakMatrix",groupBy = "Celltype",bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1",transpose = TRUE,plotLog2FC = TRUE)
p <- plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj,addDOC = FALSE)
pdf(file="Heatmap_peak_celltype.pdf", width = 8, height = 6)
p
dev.off()

#RXFP2+FIB vs FIB 
markerTest1 <- getMarkerFeatures(ArchRProj = proj,useMatrix = "PeakMatrix",groupBy ="Celltype",testMethod = "wilcoxon",bias = c("TSSEnrichment", "log10(nFrags)"),useGroups = "RXFP2+_FIB",bgdGroups = "FIB")
pma <- markerPlot(seMarker = markerTest1, name = "RXFP2+FIB_vs_FIB", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
pdf(file="MA_peak.pdf", width = 8, height = 6)
pma
dev.off()

#motif annotation
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "JASPAR2020", name = "Motif",force = TRUE)
enrichMotifs <- peakAnnoEnrichment(seMarker = markerTest1, cutOff = "FDR <= 0.05 & Log2FC >= 1",ArchRProj = proj,peakAnnotation = "Motif")
df <- data.frame(TF= rownames(enrichMotifs),mlog10Padj = assay(enrichMotifs))
df <- df[order(df$RXFP2+_FIB, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, RXFP2+_FIB, color = RXFP2+_FIB))+geom_point(size = 1)+ggrepel::geom_label_repel(data = df[rev(seq_len(20)), ], aes(x = rank, y = RXFP2+_FIB, label = TF),size = 1.5,nudge_x = 2,color = "black",max.overlaps=40)+ theme_ArchR() +ylab("-log10(P-adj) Motif Enrichment")
pdf(file="enrichMotifs.pdf", width = 4, height = 5)
ggUp
dev.off()

#save project
saveArchRProject(ArchRProj = proj,outputDirectory = './save_proj',load = FALSE,overwrite = TRUE)













































