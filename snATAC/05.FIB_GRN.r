#1.combine snRNA-snATAC
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
source("./utils.R")
source("./optMatching_functions.R") 
coembed <- readRDS("./coembed.annotation.Rds")
coembed <- subset(coembed, annotation %in% c("FIB1", "FIB2", "FIB6"))
coembed.sub <- coembed
obj.atac <- subset(coembed.sub, tech == "ATAC")
obj.rna <- subset(coembed.sub, tech == "RNA")
cca_umap_df <- as.data.frame(coembed.sub@reductions$umap_harmony_v3@cell.embeddings)
colnames(cca_umap_df) <- c("UMAP1", "UMAP2")
options(repr.plot.height = 6, repr.plot.width = 6)

df_cell_pairing <- cell_pairing(ATACpcs = obj.atac@reductions$harmony@cell.embeddings,RNApcs= obj.rna@reductions$harmony@cell.embeddings,cca_umap_df = cca_umap_df, nCores = 80)

sel_cells <- c(df_cell_pairing$ATAC,df_cell_pairing$RNA)
coembed.sub <- coembed[,sel_cells]
options(repr.plot.height = 5, repr.plot.width = 10)
p <- DimPlot(coembed.sub, reduction = "umap_harmony_v2", group.by = "annotation", split.by = "tech")
pdf("umap.pdf")
p
dev.off()
df_cell_pairing$cell_name <- paste0("cell-", 1:nrow(df_cell_pairing))

saveRDS(df_cell_pairing, "/storage/public/home/2022050408/scATAC/ATAC_RNA_matching.rds")



#2.trajectory
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(viridis))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(destiny))
suppressMessages(library(plotly))
coembed <- readRDS("./coembed.annotation.Rds")
unique(coembed$annotation)
df <- coembed@meta.data %>%
    as.data.frame() %>%
    subset(., annotation %in% c("FIB1", "FIB2", "FIB6"))
coembed <- coembed[, rownames(df)]
source("./trajectory_ArchR.R")
source("visualization.R")
source("getQuantiles.R")
traj1 <- c("FIB1", "FIB2", "FIB6")
coembed <- RunUMAP(coembed,
               dims = 1:50,
               reduction = 'harmony',
               reduction.name = "umap_harmony_v4",
               reduction.ke = 'umap_harmony_v4_',
              verbose = FALSE,
                    min.dist = 0.4)
matDR <- coembed@reductions$umap_harmony_v4@cell.embeddings
dm <- DiffusionMap(as.matrix(matDR),
                   verbose = TRUE)
as.data.frame(dm)[, c("DC1", "DC2")]
embedding <- as.data.frame(dm)[, c("DC1", "DC2")]
colnames(embedding) <- paste0("DC_", 1:ncol(embedding))
coembed[['dm']] <- CreateDimReducObject(embeddings = as.matrix(embedding),
                                              key = "DC")

cols <- ArchR::paletteDiscrete(unique(coembed@meta.data[, "annotation"]))

coembed
table(coembed$annotation)

p1 <- DimPlot(coembed, reduction = "dm", cols = cols) +xlab("DC 1") + ylab("DC 2") +theme(axis.ticks = element_blank(),axis.text = element_blank())
p1

p2 <- TrajectoryPlot(coembed, trajectory = "mesoderm_derived1", reduction = "dm", size = 1, addArrow =FALSE) +xlab("DC 1") + ylab("DC 2") +theme_cowplot() + ggtitle("mesoderm_derived")
p2


pdf("umap_trajectory_dm.pdf")
p1
dev.off()
p2
pdf("umap_trajectory_dm_fib126.pdf")
p2
dev.off()
saveRDS(coembed, file = "./coembed.annotation.trajectory.Rds")




#3.subset_ATAC
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(parallel))
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)

df.matching <- readRDS("./ATAC_RNA_matching.rds")
coembed <- readRDS("./coembed.annotation.trajectory.Rds")
obj.atac <- coembed[, df.matching$ATAC]
obj.atac <- coembed[, df.matching$mesoderm_derived1]
proj <- loadArchRProject("../subfib/", showLogo = FALSE)
proj <- subsetArchRProject(proj,
                           cells = colnames(obj.atac),
                           outputDirectory = "../trajectory_mesoderm_derived1",
                           force = TRUE)
meta.data <- as.data.frame(obj.atac@meta.data)
meta.data <- meta.data[proj@cellColData@rownames, ]
annotation <- meta.data$annotation
proj <- addCellColData(proj,
                       data = as.character(annotation),
                        cells = rownames(meta.data),
                       name = "annotation",
                       force = TRUE)
trajectory <- meta.data$mesoderm_derived1
proj <- addCellColData(proj,
                       data = as.numeric(trajectory),
                        cells = rownames(meta.data),
                       name = "mesoderm_derived1",
                       force = TRUE)
					   
embedding <- obj.atac@reductions$harmony@cell.embeddings
embedding <- embedding[rownames(proj), ]
proj@reducedDims[["Harmony"]] <- SimpleList(matDR = as.data.frame(embedding),
                                      params = NULL,
                                           date = Sys.time(),
    scaleDims = NA, #Do not scale dims after
    corToDepth = NA)
embedding <- obj.atac@reductions$umap_harmony_v4@cell.embeddings
embedding <- embedding[rownames(proj), ]
colnames(embedding) <- c("Harmony#UMAP_Dimension_1",
                         "Harmony#UMAP_Dimension_2")
proj@embeddings[["umap"]] <- SimpleList(df = as.data.frame(embedding),
                                      params = NULL)

embedding <- obj.atac@reductions$dm@cell.embeddings
embedding <- embedding[rownames(proj), ]
colnames(embedding) <- c("Harmony#DM_Dimension_1",
                         "Harmony#DM_Dimension_2")

proj@embeddings[["dm"]] <- SimpleList(df = as.data.frame(embedding),
                                      params = NULL)


pathToMacs2 <- findMacs2()
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "annotation", force = TRUE)
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "annotation",
    pathToMacs2 = pathToMacs2,genomeSize = 2.7e+09
)
proj <- addPeakMatrix(proj, binarize = TRUE, force = TRUE)
proj1 <- addMotifAnnotations(ArchRProj = proj, motifSet = "encode", name = "Motif",species = getGenome(ArchRProj = proj), force = TRUE)
proj1 <- addBgdPeaks(proj1, force = TRUE)
proj1 <- addDeviationsMatrix(
  ArchRProj = proj1,
  peakAnnotation = "Motif",
  force = TRUE,
    binarize = TRUE
)
df.matching.sub <- df.matching %>%
    dplyr::filter(ATAC %in% colnames(obj.atac))
obj.rna <- subset(coembed[, df.matching.sub$RNA])
obj.rna
geneMatrix <- getMatrixFromProject(ArchRProj = proj1,
                                   useMatrix = "GeneScoreMatrix")
gex.mat <- as.matrix(obj.rna@assays$RNA@counts)
colnames(gex.mat) <- df.matching.sub$ATAC

rowRanges <- rowData(geneMatrix)
sel_genes <- intersect(rownames(gex.mat), rowRanges$name)

length(sel_genes)

gex.mat <- gex.mat[sel_genes, ]
rownames(rowRanges) <- rowRanges$name
rowRanges <- rowRanges[sel_genes, ]

rowRanges <- GRanges(rowRanges$seqnames,
                     IRanges(start = as.numeric(rowRanges$start),
                             end = as.numeric(rowRanges$start) + 1))

seRNA <- SummarizedExperiment(assays = SimpleList(counts = gex.mat),
                              rowRanges = rowRanges)

proj1 <- addGeneExpressionMatrix(proj1,
                                seRNA = seRNA,
                                force = TRUE)
								
saveRDS(obj.atac, file = "./snATAC.trajectory.Rds")
saveArchRProject(ArchRProj = proj1,load = FALSE)


#4. tf analysis
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(parallel))
suppressMessages(library(ggrepel))
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)

p <- plotTrajectory(proj, trajectory = "mesoderm_derived1",
                    colorBy = "cellColData",
                    continuousSet = "blueYellow",
                    name = "mesoderm_derived1",
                    embedding = "dm",
                    plotAs = "points",
                    rastr = FALSE,
                    addArrow = FALSE,
                    size = 0.5)
p1 <- p[[1]] +
    theme_cowplot() +
    xlab("DC 1") + ylab("DC 2") +
    ggtitle("")

p2 <- plotEmbedding(proj, colorBy = "cellColData",
                   name = "annotation",
                    embedding = "dm",
                   plotAs = "points",
                   rastr = FALSE,
                   size = 1) + theme_cowplot() +
    xlab("DC 1") + ylab("DC 2") +
    ggtitle("")
options(repr.plot.width = 10, repr.plot.height = 5)

p1 + p2
pdf("dm_umap_trajectroy.pdf")
p1+p2
dev.off()

trajMM  <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "MotifMatrix",
                         log2Norm = FALSE,
                         scaleTo = NULL,
                        smoothWindow = 11)
trajMM <- trajMM[!grepl("deviations", rownames(trajMM)), ]
p1 <- plotTrajectoryHeatmap(trajMM,
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2))
p1

trajGEX <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "GeneExpressionMatrix",
                         log2Norm = TRUE,
                        smoothWindow = 11)
p2 <- plotTrajectoryHeatmap(trajGEX,
                        varCutOff = 0.5,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2))

p2

source("correlation.R")
df1 <- correlation_analysis(trajMM, trajGEX, cor_threshold = 0.1, varCutOff1 = 0, varCutOff2 =0.5)
0.5)
write.csv(df1, "./sel_tf_by_expression.csv")
saveArchRProject(ArchRProj = proj, load = FALSE)	



#5. peak to gene
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(parallel))
suppressMessages(library(ggrepel))
addArchRThreads(threads = 4)
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)

proj <- loadArchRProject("./", showLogo = FALSE)	

library(Rcpp)
Rcpp::sourceCpp(code='
  #include <Rcpp.h>

  using namespace Rcpp;
  using namespace std;

  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }

    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }

    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }

    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();

    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }

    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }

    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }

    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }

    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);    } 

    return(cor);

  }'
)


Rcpp::sourceCpp(code='
  #include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::IntegerVector determineOverlapCpp(IntegerMatrix m, int overlapCut){

  int k2 = 2 * m.ncol();
  int nr = m.nrow();
  int nUnion;
  int maxOverlap;
  IntegerVector unionVector;
  IntegerVector testVector = IntegerVector(nr);
  IntegerVector nOverlap = IntegerVector(nr);
  NumericVector maxOverlapVector = NumericVector(nr);
  IntegerVector vi;
  IntegerVector vj;

  for (int i = 1; i < nr; i++){
   
    if (i % 500 == 0) Rcpp::Rcout << "Completed Computing KNN Overlap " << i << " of " << nr << endl;
    
    for(int j = 0; j < i; j++){
      
      if(testVector(j) == 0){
        vi = m(i, _);
        vj = m(j, _);
        unionVector = union_( vi , vj );
        nUnion = unionVector.size();
        nOverlap(j) = k2 - nUnion;
      }else{
        nOverlap(j) = 0;
      }
    }

    maxOverlap = max( nOverlap );
    maxOverlapVector(i) = maxOverlap;
    if(maxOverlap > overlapCut){
      testVector(i) = -1;
    }

  }

  return testVector;

}'
)

trajGEX <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "GeneExpressionMatrix",
                         log2Norm = TRUE,
                        smoothWindow = 11)

ht <- plotTrajectoryHeatmap(trajGEX,
                        varCutOff = 0.9,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2))

options(repr.plot.height = 6, repr.plot.width = 6)
saveRDS(trajGEX, "./trajGEX.Rds")
ht

geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
groupMatRNA <- plotTrajectoryHeatmap(trajGEX,
                        varCutOff = 0.9,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2),
geneSet <- geneMatrix@elementMetadata
rownames(geneSet) <- paste0(geneSet$seqnames, ":", geneSet$name)
geneSet <- geneSet[rownames(groupMatRNA), ]
geneStart <- GRanges(geneSet$seqnames,IRanges(geneSet$start,width = 1),name = geneSet$name,idx= geneSet$idx)
seRNA <- SummarizedExperiment(assays = SimpleList(RNA = groupMatRNA),
                              rowRanges = geneStart)
seRNA

trajPM <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "PeakMatrix",
                         log2Norm = TRUE,
                        smoothWindow = 11)

trajPM

groupMatATAC <- plotTrajectoryHeatmap(trajPM,
                        varCutOff = 0,
                                      maxFeatures = 197282,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2),
                           returnMatrix = TRUE
                          )
nrow(groupMatATAC)
head(groupMatATAC)		

df_peak <- stringr::str_split_fixed(rownames(groupMatATAC), ":|-|_", 3)
head(df_peak)

peakSet <- GRanges(df_peak[, 1], IRanges(start = as.numeric(df_peak[, 2]),
                                           end = as.numeric(df_peak[, 3])))

seATAC <- SummarizedExperiment(assays = SimpleList(ATAC = groupMatATAC),
                               rowRanges = peakSet)
seATAC
seRNA
maxDist = 250000

o <- DataFrame(findOverlaps(resize(seRNA, 2 * maxDist + 1, "center"),
                            resize(rowRanges(seATAC), 1, "center"),
                            ignore.strand = TRUE))

o$distance <- distance(rowRanges(seRNA)[o[, 1]],
                       rowRanges(seATAC)[o[, 2]])
colnames(o) <- c("B", "A", "distance")
df <- rowRanges(seATAC)[o$A, ]
o$gene <- rowData(seRNA)[o$B, ]$name
o$peak <- paste0(df@seqnames, "_",
                 as.data.frame(df@ranges)$start, "_",
                 as.data.frame(df@ranges)$end)
outATAC <- file.path(getOutputDirectory(proj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
outRNA <- file.path(getOutputDirectory(proj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")

o$Correlation <- rowCorCpp(as.integer(o$A),as.integer(o$B),assay(seATAC),assay(seRNA))
o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) -	2)
o$FDR <- p.adjust(o$Pval, method = "fdr")
out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA",
        "VarAssayB", "distance")]
colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR",
        "VarQATAC", "VarQRNA", "Distance")
mcols(peakSet) <- NULL
names(peakSet) <- NULL
metadata(out)$peakSet <- peakSet
metadata(out)$geneSet <- geneStart
out$gene <- o$gene
out$peak <- o$peak
out <- out[!is.na(out$FDR), ]
metadata(out)$seATAC <- outATAC
metadata(out)$seRNA <- outRNA
metadata(proj@peakSet)$Peak2GeneLinks <- out
saveRDS(seATAC, file = outATAC)
saveRDS(seRNA, file = outRNA)
df_p2g <- out %>%
    as.data.frame() %>%
    subset(Correlation > 0) %>%
    subset(Distance > 2000 & FDR < 1e-04)
df_peak <- stringr::str_split_fixed(df_p2g$peak, "_", 3)
df_p2g$peak_name <- paste0(df_peak[, 1], ":", df_peak[, 2], "_", df_peak[, 3])
df_p2g$gene_name <- paste0(df_peak[, 1], ":", df_p2g$gene)
trajPM2 <- trajPM[df_p2g$peak_name, ]
trajGEX2 <- trajGEX[df_p2g$gene_name, ]
saveRDS(trajPM2, "./p2g_trajPM2.Rds")
saveRDS(trajGEX2, "./p2g_trajGEX2.Rds")		
ht1 <- plotTrajectoryHeatmap(trajPM2,
                            pal = paletteContinuous(set = "blueYellow"),
                            varCutOff = 0,
                            #limits = c(-2, 2),
                             maxFeatures = 100000,
                             rowOrder = 1:nrow(df_p2g))

ht2 <- plotTrajectoryHeatmap(trajGEX2,
                            pal = paletteContinuous(set = "horizonExtra"),
                            varCutOff = 0,
                            limits = c(-2, 2),
                             maxFeatures = 100000,
                             rowOrder = 1:nrow(df_p2g))

options(repr.plot.height = 6, repr.plot.width = 10)

ht1 + ht2
write.csv(df_p2g, file = "./p2g.csv")
saveArchRProject(ArchRProj = proj,
                 load = FALSE)
				 



#6. construct GRN
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(parallel))
suppressMessages(library(ggrepel))
suppressMessages(library(circlize))
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)

proj <- loadArchRProject("./", showLogo = FALSE)
trajMM  <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "MotifMatrix",
                         log2Norm = FALSE,
                         scaleTo = NULL,
                        smoothWindow = 11)
trajMM <- trajMM[!grepl("deviations", rownames(trajMM)), ]

p1 <- plotTrajectoryHeatmap(trajMM,
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2))

options(repr.plot.height = 8, repr.plot.width = 8)

p1
matMM <- plotTrajectoryHeatmap(trajMM,
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                              returnMatrix = TRUE)
nrow(matMM)

head(matMM)
df_tf_time_point <- data.frame(TF = rownames(matMM),
                              time_point = seq(1, 100, length.out = 2065))
rownames(df_tf_time_point) <- rownames(matMM)
head(df_tf_time_point)
df_tf <- read.csv("./sel_tf_by_expression.csv",
                  row.names = 1)
				  
df_tf_time_point <- df_tf_time_point[df_tf$name1, ]
df_tf$time_point <- df_tf_time_point$time_point
df_tf <- df_tf[order(df_tf$time_point), ]
head(df_tf)

df_tf_sub <- df_tf[, c("matchname1", "time_point")]
colnames(df_tf_sub) <- c("tf", "time_point")

p <- ggplot(data = df_tf_sub, aes(x = reorder(tf, time_point), y = time_point, label = tf)) +
    geom_point() +
    geom_text_repel() +
    xlab("TF") + ylab("Pesuto-time") +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank())
p
write.csv(df_tf_sub, "./time_point.csv")

trajMM  <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "MotifMatrix",
                         log2Norm = FALSE,
                         scaleTo = NULL,
                        smoothWindow = 11)

trajMM <- trajMM[!grepl("deviations", rownames(trajMM)), ]

trajGEX <- getTrajectory(ArchRProj = proj,
                         name = "mesoderm_derived1",
                         useMatrix = "GeneExpressionMatrix",
                         log2Norm = TRUE,
                        smoothWindow = 11)

tf_activity <- plotTrajectoryHeatmap(trajMM,
                            varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2),
                           returnMatrix = TRUE)

gene_expression <- plotTrajectoryHeatmap(trajGEX,
                        varCutOff = 0.9,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2),
                                         returnMatrix = TRUE)
head(tf_activity)
head(gene_expression)
tf_activity <- tf_activity[df_tf$name1, ]‘
df_p2g <- read.csv("./p2g.csv", row.name = 1)
sel_genes <- intersect(rownames(gene_expression), unique(df_p2g$gene_name))
length(sel_genes)
gene_expression <- gene_expression[sel_genes, ]
df_cor <- t(cor(t(tf_activity), t(gene_expression)))
df_cor <- df_cor[, df_tf$name1]
col_fun = colorRamp2(df_tf$time_point,
                     ArchR::paletteContinuous(set = "blueYellow", n = length(df_tf$time_point)))

column_ha <- HeatmapAnnotation(time_point = df_tf$time_point,
                              col = list(time_point = col_fun))

ht1 <- Heatmap(as.matrix(df_cor),
               name = "correlation",
               cluster_columns = FALSE,
             clustering_method_rows = "ward.D2",
               top_annotation = column_ha,
              show_row_names = FALSE,
              row_km = 3,
              border = TRUE)

options(repr.plot.width = 12, repr.plot.height = 10)

draw(ht1)
gene_cluster <- row_order(ht1)

gene_c1 <- rownames(df_cor[gene_cluster[[1]], ])
gene_c2 <- rownames(df_cor[gene_cluster[[2]], ])
gene_c3 <- rownames(df_cor[gene_cluster[[3]], ])

head(gene_c1)
gene_c3
length(gene_c1)
length(gene_c2)
length(gene_c3)
df_gene <- lapply(1:length(gene_cluster), function(i){
    df <- rownames(df_cor[gene_cluster[[i]], ]) %>%
        as.data.frame()

    colnames(df) <- "gene"
    df$cluster <- i

    return(df)

}) %>% Reduce(rbind, .)

df_gene$cluster <- stringr::str_replace_all(df_gene$cluster, c("1" = "FIB1",
                                                              "2" = "FIB2",
                                                              "3" = "FIB6"))
df_gene$gene <- stringr::str_split_fixed(df_gene$gene, ":", 2)[, 2]
write.csv(df_gene, "./gene_cluster.csv")		

mat_p2g <- df_p2g %>%
    select(c(peak_name, gene_name, Correlation)) %>%
    tidyr::pivot_wider(names_from = peak_name, values_from = Correlation) %>%
    textshape::column_to_rownames("gene_name")
mat_p2g <- df_p2g %>%
    dplyr::select(c(peak_name, gene_name, Correlation)) %>%
    tidyr::pivot_wider(names_from = peak_name, values_from = Correlation) %>%
    textshape::column_to_rownames("gene_name")
mat_p2g[is.na(mat_p2g)] <- 0
mat_p2g[mat_p2g>0] <- 1

motif_matching <-getMatches(proj)
rowRanges <- rowRanges(motif_matching)
chr <- seqnames(rowRanges)
ranges <- as.data.frame(ranges(rowRanges))

matches <- as.matrix(assays(motif_matching)$matches)
colnames(matches) <- paste0("z:", colnames(matches))
rownames(matches) <- paste0(chr, ":", ranges$start, "_", ranges$end)
matches <- matches[colnames(mat_p2g),
                   colnames(df_cor)]
head(matches)

gene_tf <- as.matrix(mat_p2g) %*% as.matrix(matches)
head(gene_tf)

gene_tf[gene_tf>0] <- 1
df_cor <- df_cor * gene_tf
suppressMessages(library(igraph))
df_cor <- as.data.frame(df_cor)
df_cor$gene <- rownames(df_cor)
df_cor_2 <- df_cor %>%
    tidyr::pivot_longer(!gene, names_to = "tf", values_to = "correlation") %>%
    subset(correlation > 0.5) %>%
    dplyr::select(c(tf, gene, correlation))
df_cor_2 <- df_cor_2[!grepl("-AS", df_cor_2$gene), ]
df_cor_2$tf <- stringr::str_split_fixed(df_cor_2$tf, ":", 2)[, 2]
df_cor_2$tf <- stringr::str_split_fixed(df_cor_2$tf, "_", 2)[, 1]
df_cor_2$gene <- stringr::str_split_fixed(df_cor_2$gene, ":", 2)[, 2]

write.csv(df_cor_2, "./gene_regulatory_network.csv", 
          row.names = FALSE)		

g <- graph_from_data_frame(df_cor_2, directed=FALSE)

tf_colors <- ArchR::paletteContinuous(set = "blueYellow", n = 65)
names(tf_colors) <- df_tf$name1

gene_colors <- rep("white", length(unique(df_cor_2$gene)))
names(gene_colors) <- unique(df_cor_2$gene)

vertex.color <- c(tf_colors, gene_colors)


tf_colors <- ArchR::paletteContinuous(set = "blueYellow", n = 65)
names(tf_colors) <- df_tf$name1

gene_colors <- rep("black", length(unique(df_cor_2$gene)))
names(gene_colors) <- unique(df_cor_2$gene)

vertex.frame.color <- c(tf_colors, gene_colors)

co <- layout_with_fr(g, dim = 2, niter = 1000)

## only show labels for TF
V(g)$label <- ifelse(V(g)$name %in% df_cor_2$tf, V(g)$name, NA)

## for TFs, we increase the size of nodes
V(g)$size <- ifelse(V(g)$name %in% df_cor_2$tf, 3, 1)

E(g)$weight <- E(g)$correlation

options(repr.plot.height = 25, repr.plot.width = 25)

plot(g, layout = co, 
     vertex.color = vertex.color,
     vertex.frame.color = vertex.frame.color,
     vertex.label.dist = 0.5,
     edge.width = E(g)$weight*2,
     vertex.label.cex = 2,
     edge.color = adjustcolor("gray", alpha = 1),
     edge.curved=seq(-1, 1, length = ecount(g))
    )
sel.tfs <- intersect(df_cor_2$tf, df_cor_2$gene)
sel.tfs

df_cor_3 <- subset(df_cor_2, gene %in% sel.tfs)

g <- graph_from_data_frame(df_cor_3, directed=TRUE)
co <- layout_with_fr(g, dim = 2, niter = 1000)

options(repr.plot.height = 20, repr.plot.width = 20)

plot(g, layout = co,
    vertex.label.cex = 2,
    vertex.size = 7,
         edge.color = adjustcolor("gray", alpha = 1)
    )
write.csv(df_cor_3, "./tf_tf_network.csv", row.names = FALSE)


pagerank <- page_rank(g, weights = E(g)$weights)
bet <- betweenness(g,weights = E(g)$weights, normalized = TRUE)
df_measure <- data.frame(
    tf = V(g)$name,
    pagerank = pagerank$vector,
    betweenness = bet
  ) 
wirte.csv(df_measure,"netmeasures_final_cms.csv")








# 7.GRN plot
suppressMessages(library(ArchR))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(circlize))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidygraph))
suppressMessages(library(patchwork))
library(BSgenome.Btaurus.UCSC.bosTau9)
library(TxDb.Btaurus.UCSC.bosTau9.refGene)
library(org.Bt.eg.db)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Btaurus.UCSC.bosTau9)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Btaurus.UCSC.bosTau9.refGene, OrgDb = org.Bt.eg.db)

df <- read.csv("./gene_regulatory_network.csv")
tf.list <- unique(df$tf)
gene.list <- setdiff(unique(df$gene), tf.list)

length(tf.list)
length(gene.list)

head(tf.list)
head(gene.list)

df_measure <-  read.csv("./netmeasures_final_cms.csv", row.names = 1)  %>%
        mutate(pagerank = scale(pagerank)[, 1]) %>%
        mutate(betweenness = scale(betweenness)[, 1])
min.page <- min(df_measure$pagerank)
min.bet <- min(df_measure$betweenness)
df_measure$importance <-
    sqrt((df_measure$pagerank - min.page) ** 2 +
           (df_measure$betweenness - min.bet) ** 2)
df_measure <- df_measure[order(-df_measure$importance), ]
rownames(df_measure) <- df_measure$tf

p1 <- ggplot(data = df_measure, aes(x = pagerank, y = betweenness)) +
    geom_point() +
    xlab("centrality") + ylab("betweenness") +
    geom_text_repel(aes(label = tf)) +
    theme_cowplot()
p2 <- ggplot(data = df_measure, aes(x = reorder(tf, -importance), y = importance)) +
    geom_point() +
    xlab("Rank") + ylab("Importance") +
    geom_text_repel(aes(label = tf)) +
    geom_hline(yintercept = 2, color = "red") +
    theme_cowplot() +
    theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank())

p1 + p2

df_gene_clustering <- read.csv("./gene_cluster.csv",
                         row.names = 1) %>%
    subset(gene %in% gene.list)

df_timepoint <- read.csv("./time_point.csv",
                        row.names = 1) %>%
    subset(tf %in% tf.list)

g <- graph_from_data_frame(df, directed=TRUE)
names(tf_size) <- df_measure$tf
names(gene_size) <- gene.list
tf_size <- df_measure$importance
gene_size <- rep(min(df_measure$importance), length(gene.list))
v_size <- c(tf_size, gene_size)
V(g)$size <- v_size[V(g)$name]

df_timepoint <- df_timepoint[order(df_timepoint$time_point), ]
head(df_timepoint)
tf_color <- ArchR::paletteContinuous(set = "blueYellow",
                                     n = nrow(df_timepoint))
names(tf_color) <- df_timepoint$tf
gene_color <- stringr::str_replace_all(df_gene_clustering$cluster,
                                      c("FIB1" = "#D51F26",
                                       "FIB2" = "#208A42",
                                       "FIB6" = "#272E6A"))
names(gene_color) <- df_gene_clustering$gene
v_color <- c(tf_color, gene_color)

tf_alpha <- rep(1, length(tf.list))
gene_alpha <- rep(0.5, length(gene.list))
names(tf_alpha) <- tf.list
names(gene_alpha) <- gene.list

v_alpha <- c(tf.list, gene.list)
V(g)$alpha <- v_alpha[V(g)$name]

layout <- layout_with_fr(g, weights = E(g)$correlation, dim = 2, niter = 1000)
df_measure_sub <- subset(df_measure, importance > 0)
df <- subset(df, tf %in% df_measure_sub$tf)
write.csv(df, "./gene_regulatory_network_filtered.csv")

ggraph(g, layout = layout) +
geom_edge_link(edge_colour = "gray", edge_alpha = 0.25) +
geom_node_point(aes(size = V(g)$size,
                   color = as.factor(name),
                   alpha = V(g)$alpha,
                   ),
                show.legend = FALSE,max.overlaps = Inf) +
scale_size(range = c(1, 10)) +
scale_color_manual(values = v_color) +
geom_node_label(aes(filter = V(g)$name %in% df_measure_sub$tf,
                    label = V(g)$name),
                    repel = TRUE,
                    hjust = "inward",
                size = 5,
                show.legend = FALSE,max.overlaps = Inf) +
theme_void()


	