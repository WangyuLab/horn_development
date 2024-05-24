###  differential expressed genes(DEGs) of cell types###

# load packages
rm(list=ls())
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(paletteer)  
library(gplots)

# read the data
sce.all=readRDS( "../obj2.rds")

# divide cell types to cell lists
cell_list = split(colnames(sce.all),sce.all$celltype)
cell_list
names(cell_list)

# use function FindAllMarkers, export top10 markers dotplot and heatmap, save Rdata of top markers
for ( pro in names(cell_list) ) {
  sce=sce.all[,colnames(sce.all) %in% cell_list[[pro]]]
  sce <- CreateSeuratObject(counts = sce@assays$RNA@counts, 
                            meta.data = sce@meta.data, 
                            min.cells = 3, 
                            min.features = 200)  
  sce <- NormalizeData(sce)  
  sce = FindVariableFeatures(sce)
  sce = ScaleData(sce, 
                  vars.to.regress = c("nFeature_RNA",
                                      "percent_mito"))
  
  Idents(sce)=sce$orig.ident
  table(Idents(sce))
  sce.markers <- FindAllMarkers(object = sce, 
                                only.pos = TRUE, 
                                logfc.threshold = 0.1, 
                                min.pct = 0.1, 
                               thresh.use = 0.1)
  
  write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
  sce.markers=sce.markers[order(sce.markers$cluster,
                                sce.markers$avg_log2FC),]
  top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  
  sce.Scale <- ScaleData( sce ,features =  unique(top10$gene)  )  
  
  DoHeatmap(sce.Scale,
            features =  unique(top10$gene),
            assay = 'RNA', label = T)+
    scale_fill_gradientn(colors = c("white","grey","firebrick3"))
  
  ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),
         height = 8)
  
  p <- DotPlot(sce , features = unique(top10$gene)  ,
               assay='RNA'  )  + coord_flip()
  
  p
  ggsave(plot=p, 
         filename=paste0("check_top10-marker_by_",
                         pro,"_cluster.pdf") 
         ,height = 8)
  save(sce.markers,
       file=paste0(pro,'_sce.markers.Rdata'))	


  }	   