### construct single cell trajectories using monocle2 pseudotime analysis  ###

# load packages
library(Seurat)
library(dplyr)
library(Matrix)
library(patchwork)
library(monocle)
library(ggplot2)

# read the data
neu <- readRDS("neu.rds")

# extract phenotypic and genic information
data <- as(as.matrix(neusub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = neusub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
row.names(data)

# construct a CDS object
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
# estimate size factors and dispersions							  
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

# filter low quality cells
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >=10))

# choose genes that define a cell's progress (four methods)
deg.cluster <- FindAllMarkers(neu)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
monocle_cds <- setOrderingFilter(monocle_cds, express_genes)

# reduce data dimensionality
monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
    method = 'DDRTree')

# order cells along the trajectory
monocle_cds <- orderCells(monocle_cds)

# visualize the trajectory 
plot_cell_trajectory(monocle_cds, color_by = "neu_celltype")
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
plot_cell_trajectory(monocle_cds, color_by = "State")

# differential expression analysis of branch
BEAM_res <- BEAM(cds, branch_point = 2, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_genes <- top_n(BEAM_res,n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,], branch_point = 2, num_clusters =4,show_rownames = T,return_heatmap = T)
ggsave("BEAM_heatmap.pdf",p$ph_res, width = 6.5, height = 10)

# save image and rds
saveRDS(neu, "neu.rds")
save.image("neu_monocle2.rd")
savehistory("neu_monocle2.r")