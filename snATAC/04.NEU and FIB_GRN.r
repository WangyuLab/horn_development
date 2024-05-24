Neural cells and fibroblasts were extracted for subcluster analysis

# Use scMEGA to integrate snATAC-seq and snRNA-seq data from three periods
# D100
ATAC100 <- readRDS('{file_pathway}/scATAC-seq/horn_200d/save_proj/Save-ArchR-Project.rds')
peak_counts <- readRDS(glue::glue("{file_pathway}/scATAC-seq/horn_100d/peak_counts.Rds"))
metadata <- as.data.frame(ATAC100@cellColData) %>%subset(., select = c("Sample","Clusters"))
peak_counts_sub <- peak_counts[,colnames(peak_counts) %in% rownames(metadata)]
chrom_assay <- CreateChromatinAssay(counts = peak_counts_sub,sep = c(":", "-"),min.cells = 100)
obj.atac<-CreateSeuratObject(counts=chrom_assay,assay="ATAC",meta.data=metadata,names.field=1,names.delim="#")
embedding <- ATAC100@embeddings$UMAP$df
colnames(embedding) <- paste0("UMAP_", 1:ncol(embedding))
embedding <- embedding[colnames(obj.atac), ]
obj.atac[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding),assay = "ATAC",key ="UMAP_")
obj.rna <- readRDS("{file_pathway}/snRNA/horn/100dhorn.rds")

#
atac <- getMatrixFromProject(ArchRProj = ATAC100,useMatrix = "GeneScoreMatrix")
gene_counts <- atac@assays@data$GeneScoreMatrix
rownames(gene_counts) <- atac@elementMetadata$name
saveRDS(gene_counts, file = "gene_activity100.Rds")
gene.activity <- readRDS(glue::glue("{file_pathway}/gene_activity100.Rds"))
gene.activity <- gene.activity[, colnames(obj.atac)]

#
obj.coembed<-CoembedData(obj.rna,obj.atac, gene.activity, weight.reduction ="umap", dims = 1:2,verbose =TRUE)
obj.coembed<-RunHarmony(obj.coembed,group.by.vars =c( "tech"),reduction ="pca",max.iter.harmony =30,dims.use =1:30,project.dim =FALSE,plot_convergence =FALSE)
obj.coembed<-RunUMAP(obj.coembed,dims =1:30,reduction ='harmony',reduction.name ="umap_harmony",reduction.ke ='umapharmony_',verbose =FALSE,min.dist =0.4)
p<-DimPlot(obj.coembed, group.by ="tech", reduction ="umap_harmony")
pdf(file="UMAP100_co_har_tech.pdf", width = 8, height = 6)
p
dev.off()

# first use high resolution to get a large number of clusters
obj.coembed <- FindNeighbors(obj.coembed, reduction = "harmony", dims = 1:30)
obj.coembed <- FindClusters(obj.coembed, resolution = 1, verbose = FALSE)
cols <- ArchR::paletteDiscrete(obj.coembed@meta.data[, "RNA_snn_res.1"])
p <- CellPropPlot(obj.coembed,group.by = "tech",prop.in = "RNA_snn_res.1")
pdf(file="prop_coembed_tech.pdf", width = 8, height = 6)
p
dev.off()

#Delete clusters from extreme technical sources
Idents(obj.coembed)<-"RNA_snn_res.1"
coembed.sub<-subset(obj.coembed, idents =c(42, 49), invert =TRUE)

#Re-perform UMAP embedding and clustering at a lower resolution to reduce complexity
coembed.sub<-RunUMAP(coembed.sub, dims =1:30, reduction ='harmony',reduction.name ="umap_harmony",reduction.key ='umap_harmony_',verbose =FALSE,min.dist =0.4)
coembed.sub <- FindNeighbors(coembed.sub, reduction = "harmony", dims = 1:30)
coembed.sub <- FindClusters(coembed.sub, resolution = 0.5, verbose = FALSE)
cols <- ArchR::paletteDiscrete(coembed.sub@meta.data[, "RNA_snn_res.0.5"])
p <- DimPlot(coembed.sub, group.by = "RNA_snn_res.0.5", label = TRUE,reduction = "umap_harmony",shuffle = TRUE, cols = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="newUAMP100_select_clu.pdf", width = 8, height = 6)
p
dev.off()
p<-DimPlot(coembed.sub, group.by ="RNA_snn_res.0.5", label =TRUE,reduction ="umap_harmony", shuffle =TRUE, split.by ="tech",cols =cols)+xlab("UMAP1")+ylab("UMAP2")
pdf(file="newUMAP100_select_clu_tech.pdf", width = 12, height = 6)
p
dev.off()

#Cell pairing
df.pair <- PairCells(object = coembed.sub, reduction = "harmony",pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")
sel_cells<-c(df.pair$ATAC, df.pair$RNA)
coembed.sub2<-coembed.sub[, sel_cells]
options(repr.plot.height =5, repr.plot.width =10)
p <- DimPlot(coembed.sub2, reduction ="umap_harmony", group.by ="RNA_snn_res.0.5", split.by ="tech", cols =cols)
pdf(file="newUMAP100_pair_tech.pdf", width = 12, height = 6)
p
dev.off()

#Create a new Seurat object for those paired units as if they were produced by joint sequencing
obj.pair<-CreatePairedObject(df.pair =df.pair, object =coembed.sub2,use.assay1 ="RNA", use.assay2 ="ATAC")
saveRDS(obj.pair, "./obj.pair100.rds")

#Get 150day and 200day obj.pair using the same pipeline

#Three periods of integration of obj.pair objects
pair100 <- readRDS("{file_pathway}/obj.pair100.rds")
pair150 <- readRDS("{file_pathway}/obj.pair150.rds")
pair200 <- readRDS("{file_pathway}/obj.pair200.rds")
pair100$Sample <- "100d"
pair150$Sample <- "150d"
pair200$Sample <- "200d"
100d_150d <- merge(pair100,y =pair150,add.cell.ids = c("100d","150d"),project = "merge")
merge <- merge(merge,y =pair200,add.cell.ids = c("","200d"),project = "merge")

#Reduction and clustering
merge <- NormalizeData(merge, normalization.method = "LogNormalize", scale.factor = 10000)
merge <- FindVariableFeatures(merge,selection.method = "vst",nfeatures = 2000)
merge <- ScaleData(merge,features=VariableFeatures(merge))
merge <- RunPCA(merge, features = VariableFeatures(object = merge))
merge <- RunHarmony(merge,group.by.vars =c("Sample"),reduction ="pca",max.iter.harmony =30,dims.use =1:30,project.dim =FALSE,plot_convergence =FALSE)
merge <- RunUMAP(merge,dims =1:30,reduction ='harmony',reduction.name ="umap_harmony",reduction.ke ='umapharmony_',verbose =FALSE,min.dist =0.4)

#Next, we performed subclustering to identify different populations in our multicomponent data.
obj.coembed <- merge
coembed.sub<-RunUMAP(coembed.sub, dims =1:30, reduction ='harmony',reduction.name ="umap_harmony",reduction.key ='umap_harmony_',verbose =FALSE,min.dist =0.4)
coembed.sub <- FindNeighbors(coembed.sub, reduction = "harmony", dims = 1:30)
coembed.sub <- FindClusters(coembed.sub, resolution = 0.5, verbose = FALSE)
cols <- ArchR::paletteDiscrete(coembed.sub@meta.data[, "RNA_snn_res.0.5"])
p <- DimPlot(coembed.sub, group.by = "RNA_snn_res.0.5", label = TRUE,reduction = "umap_harmony",shuffle = TRUE, cols = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="newUAMP_select_clu0.5.pdf", width = 8, height = 6)
p
dev.off()

#Draw the first 3 markers of each cluster
all.markers <- FindAllMarkers(coembed.sub,only.pos = TRUE,min.pct = 0.25, logfc.threshold = 0.5)
df<-all.markers%>%group_by(cluster)%>%slice_max(n =3, order_by =avg_log2FC)
p<-DotPlot(coembed.sub, features =unique(df$gene))+RotatedAxis()
pdf(file="newmarker_top3marker_subclu0.5.pdf", width = 18, height = 6)
p
dev.off()
p2 <- CellPropPlot(coembed.sub,group.by = "Sample",prop.in = "RNA_snn_res.0.5")
pdf(file="newprop_clu0.5_sample.pdf", width = 8, height = 6)
p2
dev.off()

#plot select marker
p <- FeaturePlot(coembed.sub, features = c("COL1A1","COL3A1","PDGFRA","KRT15","KRT17","KRT14","PECAM1","KDR","GRIK2","PPP2R2B","ITGA7","MKI67","MYH3","TYR","MMRN1","PROX1","RXFP2"),reduction ="umap_harmony",min.cutoff = "q10", max.cutoff = "q90")
pdf(file="newUMAP_marker.pdf", width =8, height = 6)
p
dev.off()

#Identified cell type
NEU=c(13)
LEC=c(15)
END=c(3)
MEL=c(16)
MAC=c(12)
MUS1=c(10)
MUS2=c(17)
KRT=c(18,5,14,8)
FIB=c(19,20,1,7,11,6,0,2,4)
CYC=c(9)
current.cluster.ids <- c(NEU,LEC,END,MEL,MAC,MUS1,MUS2,KRT,FIB,CYC)
new.cluster.ids <- c(rep("NEU",length(NEU)),rep("LEC",length(LEC)),rep("END",,length(END)),rep("MEL",length(MEL)),rep("MAC",length(MAC)),rep("MUS1",length(MUS1)),rep("MUS2",length(MUS2)),rep("KRT",length(KRT)),rep("FIB",length(FIB)),rep("CYC",length(CYC)))
coembed.sub@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(coembed.sub@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

#plot
p<-DimPlot(coembed.sub, group.by ="celltype", reduction ="umap_harmony")
pdf(file="newUAMP_select_celltype0.5.pdf", width = 8, height = 6)
p
dev.off()
p2 <- CellPropPlot(coembed.sub,group.by = "Sample",prop.in = "celltype")+coord_flip()
pdf(file="newprop_celltype_sample.pdf", width = 6, height = 6)
p2
dev.off()

#
saveRDS(coembed.sub, "./coembed.rds")

###Take a subset of neural cells
proj <- subset(x = proj, idents = c(13))
coembed.sub<-RunUMAP(proj, dims =1:30, reduction ='harmony',reduction.name ="umap_harmony",reduction.key ='umap_harmony_',verbose =FALSE,min.dist =0.4)
coembed.sub <- FindNeighbors(coembed.sub, reduction = "harmony", dims = 1:30)
coembed.sub <- FindClusters(coembed.sub, resolution = 0.2, verbose = FALSE)
cols <- ArchR::paletteDiscrete(coembed.sub@meta.data[, "RNA_snn_res.0.2"])
p <- DimPlot(coembed.sub, group.by = "RNA_snn_res.0.2", label = TRUE,reduction = "umap_harmony",shuffle = TRUE, cols = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="UAMP_NEU_clu.pdf", width = 8, height = 6)
p
dev.off()

#show marker
all.markers <- FindAllMarkers(coembed.sub,only.pos = TRUE,min.pct = 0.5, logfc.threshold = 0.5)
df<-all.markers%>%group_by(cluster)%>%slice_max(n =3, order_by =avg_log2FC)
p<-DotPlot(coembed.sub, features =unique(df$gene))+RotatedAxis()
pdf(file="marker_top3.pdf", width = 9, height = 5)
p
dev.off()
p<-DotPlot(coembed.sub, features =c("MKI67","NGFR","SOX10","SNCA","RELN","PAX3","POSTN","RUNX2","RXFP2"),group.by = "cluster")+RotatedAxis()
pdf(file="Dotplot_clu_RNAmarker.pdf", width =9, height = 6)
p
dev.off()

# Identified cell type
Neural_crest=c(3)
sub_NEU2=c(0,2)
sub_NEU1=c(1)
current.cluster.ids <- c(Neural_crest,sub_NEU2,sub_NEU1)
new.cluster.ids <- c(rep("Neural_crest",length(Neural_crest)),rep("sub_NEU2",length(sub_NEU2)),rep("sub_NEU1",length(sub_NEU1))) coembed.sub@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(coembed.sub@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)

# plot
p2<-DimPlot(coembed.sub, group.by ="celltype", reduction ="umap_harmony")
pdf(file="UMAP_NEU_celltype.pdf", width = 8, height = 6)
p2
dev.off()

# Trajectory1 GRN
obj<-AddTrajectory(object =coembed.sub, trajectory =c(3,1),name="traj1",group.by ="RNA_snn_res.0.2", reduction ="umap_harmony",dims =1:2, use.all =TRUE)
obj <- coembed.sub[, !is.na(coembed.sub$traj1)]
p1 <- DimPlot(obj, reduction = "umap_harmony",group.by = "RNA_snn_res.0.2", cols = cols) +xlab("UMAP 1")+ ylab("UMAP 2") +ggtitle("Clusters")
p2 <- TrajectoryPlot(object = obj,trajectory = "traj1",reduction = "umap_harmony",continuousSet ="blueYellow",size = 1,addArrow = FALSE)
pdf(file="traj_1.pdf", width =10, height = 6)
p1+p2
dev.off()

# Find motif
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(object = obj,genome = BSgenome.Btaurus.UCSC.bosTau9,pfm =pfm,assay ="ATAC")
obj<-RunChromVAR(object =obj,genome =BSgenome.Btaurus.UCSC.bosTau9,assay ="ATAC")
res<-SelectTFs(object =obj, trajectory.name = "traj1",return.heatmap =TRUE,cor.cutoff =0.1,p.cutoff = 0.05)                           
df.cor<-res$tfs
ht<-res$heatmap
pdf(file="traj1_heatmap_TF2gene.pdf", width =7, height = 5)
ht
dev.off()

#
res <- SelectGenes(object = obj,trajectory.name = "traj1",var.cutoff.gene=0.75,fdr.cutoff = 1e-02,labelTop1 = 0,labelTop2 = 0)
df.p2g <- res$p2g
ht <- res$heatmap
pdf(file="traj1_heatmap_peak2gene.pdf", width =8, height = 6)
ht
dev.off()
tf.gene.cor <- GetTFGeneCorrelation(object = obj,tf.use = df.cor$tfs,gene.use = unique(df.p2g$gene),tf.assay = "chromvar",gene.assay = "RNA",trajectory.name = "traj1")
ht <- GRNHeatmap(tf.gene.cor,tf.timepoint = df.cor$time_point)
pdf(file="GRN_heatmap_traj1.pdf", width =8, height = 3)
ht
dev.off()

#
motif.matching <- obj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names
motif.matching <-motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
df.grn <- GetGRN(motif.matching = motif.matching,df.cor = tf.gene.cor,df.p2g = df.p2g)

#Visual network
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs
df.grn2 <- df.grn %>%subset(correlation > 0) %>%select(c(tf, gene, correlation)) %>%rename(weights= correlation)
p <- GRNPlot(df.grn2,tfs.timepoint = tfs.timepoint,show.tf.labels = TRUE,seed = 42,plot.importance =FALSE,min.importance = 2,remove.isolated = FALSE)
options(repr.plot.height = 20, repr.plot.width = 20)
pdf(file="GRN_traj1.pdf", width =16, height = 6)
p
dev.off()

# Trajectory2 GRN
obj<-AddTrajectory(object =coembed.sub, trajectory =c(3,0,2),name="traj2",group.by ="RNA_snn_res.0.2", reduction ="umap_harmony",dims =1:2, use.all =TRUE)
obj <- coembed.sub[, !is.na(coembed.sub$traj2)]
p1 <- DimPlot(obj, reduction = "umap_harmony",group.by = "RNA_snn_res.0.2", cols = cols) +xlab("UMAP 1")+ ylab("UMAP 2") +ggtitle("Clusters")
p2 <- TrajectoryPlot(object = obj,trajectory = "traj2",reduction = "umap_harmony",continuousSet ="blueYellow",size = 1,addArrow = FALSE)
pdf(file="traj_2.pdf", width =10, height = 6)
p1+p2
dev.off()

# Find motif
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(object = obj,genome = BSgenome.Btaurus.UCSC.bosTau9,pfm =pfm,assay ="ATAC")
obj<-RunChromVAR(object =obj,genome =BSgenome.Btaurus.UCSC.bosTau9,assay ="ATAC")
res<-SelectTFs(object =obj, trajectory.name = "traj2",return.heatmap =TRUE,cor.cutoff =0.1,p.cutoff = 0.05)                           
df.cor<-res$tfs
ht<-res$heatmap
pdf(file="traj2_heatmap_TF2gene.pdf", width =7, height = 5)
ht
dev.off()

#
res <- SelectGenes(object = obj,trajectory.name = "traj2",var.cutoff.gene=0.75,fdr.cutoff = 1e-02,labelTop1 = 0,labelTop2 = 0)
df.p2g <- res$p2g
ht <- res$heatmap
pdf(file="traj2_heatmap_peak2gene.pdf", width =8, height = 6)
ht
dev.off()
tf.gene.cor <- GetTFGeneCorrelation(object = obj,tf.use = df.cor$tfs,gene.use = unique(df.p2g$gene),tf.assay = "chromvar",gene.assay = "RNA",trajectory.name = "traj2")
ht <- GRNHeatmap(tf.gene.cor,tf.timepoint = df.cor$time_point)
pdf(file="GRN_heatmap_traj2.pdf", width =8, height = 3)
ht
dev.off()

#
motif.matching <- obj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names
motif.matching <-motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
df.grn <- GetGRN(motif.matching = motif.matching,df.cor = tf.gene.cor,df.p2g = df.p2g)

#Visual traj2 network
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs
df.grn2 <- df.grn %>%subset(correlation > 0) %>%select(c(tf, gene, correlation)) %>%rename(weights= correlation)
p <- GRNPlot(df.grn2,tfs.timepoint = tfs.timepoint,show.tf.labels = TRUE,seed = 42,plot.importance =FALSE,min.importance = 2,remove.isolated = FALSE)
options(repr.plot.height = 20, repr.plot.width = 20)
pdf(file="GRN_traj2.pdf", width =16, height = 6)
p
dev.off()
#
saveRDS(obj, "./NEU_obj.traj.rds")

###Take a subset of FIB cells
proj <- subset(x = coembed.sub, idents = c(19,20,1,7,11,6,0,2,4))

#
coembed.sub<-RunUMAP(proj, dims =1:30, reduction ='harmony',reduction.name ="umap_harmony",reduction.key ='umap_harmony_',verbose =FALSE,min.dist =0.4)
coembed.sub <- FindNeighbors(coembed.sub, reduction = "harmony", dims = 1:30)
coembed.sub <- FindClusters(coembed.sub, resolution = 0.2, verbose = FALSE)
cols <- ArchR::paletteDiscrete(coembed.sub@meta.data[, "RNA_snn_res.0.2"])
p <- DimPlot(coembed.sub, group.by = "RNA_snn_res.0.2", label = TRUE,reduction = "umap_harmony",shuffle = TRUE, cols = cols) +xlab("UMAP1") + ylab("UMAP2")
pdf(file="UAMP_FIB_clu.pdf", width = 8, height = 6)
p
dev.off()

#Draw the selected mrker
p<-DotPlot(coembed.sub, features =c("MKI67","NGFR","SOX10","SNCA","RELN","PAX3","POSTN","RUNX2","RXFP2"),group.by = "celltype")+RotatedAxis()
pdf(file="Dotplot_clu_RNAmarker.pdf", width =9, height = 6)
p
dev.off()

#draw trajectory
obj <- AddTrajectory(object =coembed.sub, trajectory =c(4,0,1,5),name="traj1",group.by ="RNA_snn_res.0.2", reduction ="umap_harmony",dims =1:2, use.all =TRUE)
obj <- obj[, !is.na(obj$traj1)]
p1 <- DimPlot(obj, reduction = "umap_harmony",group.by = "RNA_snn_res.0.2", cols = cols) +xlab("UMAP 1")+ ylab("UMAP 2") +ggtitle("Clusters")
p2 <- TrajectoryPlot(object = obj,trajectory = "traj1",reduction = "umap_harmony",continuousSet ="blueYellow",size = 1,addArrow = FALSE)
pdf(file="traj_1.pdf", width =10, height = 6)
p1+p2
dev.off()

#
saveRDS(obj, "./FIB_obj.traj.rds")

#
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
obj <- AddMotifs(object = obj,genome = BSgenome.Btaurus.UCSC.bosTau9,pfm =pfm,assay ="ATAC")
obj<-RunChromVAR(object =obj,genome =BSgenome.Btaurus.UCSC.bosTau9,assay ="ATAC")
res<-SelectTFs(object =obj, trajectory.name = "traj1",return.heatmap =TRUE,cor.cutoff =0.1)
df.cor<-res$tfs
ht<-res$heatmap
pdf(file="traj1_heatmap_TF2gene.pdf", width =7, height = 5)
ht
dev.off()

#
res <- SelectGenes(object = obj,trajectory.name = "traj1",labelTop1 = 0,labelTop2 = 0)
df.p2g <- res$p2g
ht <- res$heatmap
pdf(file="traj1_heatmap_peak2gene.pdf", width =8, height = 6)
ht
dev.off()
tf.gene.cor <- GetTFGeneCorrelation(object = obj,tf.use = df.cor$tfs,gene.use = unique(df.p2g$gene),tf.assay = "chromvar",gene.assay = "RNA",trajectory.name = "traj1")
ht <- GRNHeatmap(tf.gene.cor,tf.timepoint = df.cor$time_point)
pdf(file="GRN_heatmap_traj1.pdf", width =7, height = 3)
ht
dev.off()

#
motif.matching <- obj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names
motif.matching <-motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
df.grn <- GetGRN(motif.matching = motif.matching,df.cor = tf.gene.cor,df.p2g = df.p2g)

#
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs
df.grn2 <- df.grn %>%subset(correlation > 0.85) %>%select(c(tf, gene, correlation)) %>%rename(weights = correlation)
p <- GRNPlot(df.grn2,tfs.timepoint = tfs.timepoint,show.tf.labels = TRUE,seed = 42,plot.importance = FALSE,min.importance = 2,remove.isolated = FALSE)
options(repr.plot.height = 20, repr.plot.width = 20)
pdf(file="GRN_traj1.pdf", width =16, height = 6)
p
dev.off()































































































