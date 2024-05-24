### spatial transcriptomics data analysis using CARD ###
### cell type deconvolution ###

# install packages
install.packages('devtools')
devtools::install_github('xuranw/MuSiC')
devtools::install_github('YingMa0107/CARD')

# load packages
library(Seurat)
library(dplyr)
library(Matrix)
library(tibble)
library(ggplot2)
library(CARD)
library(MuSiC)
library(gtools)
library(scatterpie)
library(RColorBrewer)

# read the scRNA data
all <- readRDS("/home/../merge/obj5.rds")
ls <- GetAssayData(object = all, slot = "counts")
test_meta <-  data.frame(Cell = rownames(all@meta.data),cell_type = Idents(all))
head(test_meta)
test_meta$sampleInfo= "sample1"

# read the spatial transcriptomics data, including location data and coun data
loc <- read.table("/home/../visium/report/05.AllheStat/heAuto_level_matrix/subdata/L13_heAuto/barcodes_pos.tsv.gz")
names(loc) <- c("","x","y")
row.names(loc) <- loc[,1]
loc=loc[,-1]
p = ggplot(loc,aes(x=x,y=y))+geom_point(colour="#990000",size=0.1)
p + geom_vline(xintercept = c(480,770)) +geom_hline(yintercept = 470)
loc1 <- filter(loc,y>470)  ## slecet the interesting loacation
ggplot(loc1,aes(x=x,y=y))+geom_point(colour="#990000",size=0.1)

expr=Read10X("/home/yuwang/yxx/visium/report/05.AllheStat/heAuto_level_matrix/subdata/L13_heAuto/",cell.column = 1)
expr1 <- expr[,rownames(loc1)]
expr1 [1:4,1:4]

# create an CARD object
CARD_obj = createCARDObject(sc_count = ls,sc_meta = test_meta,spatial_count = expr1,spatial_location = loc1,ct.varname = "cell_type",ct.select = unique(test_meta$cell_type),sample.varname= "sampleInfo",minCountGene = 100,minCountSpot = 5)

# deconvolution using CARD
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

# visualize the proportion for each cell type
p1 <- CARD.visualize.pie(proportion = CARD_obj@Proportion_CARD,spatial_location = CARD_obj@spatial_location, colors = brewer.pal(name="Paired",11))
pdf("celltype.pdf")
p1
dev.off()

ct.visualize = c("FIB","KRT","END","MUS1","CYC","MUS2","NEU","FIB","MAC","LEC","MEL")
p2 <- CARD.visualize.prop(
	proportion = CARD_obj@Proportion_CARD,        
	spatial_location = CARD_obj@spatial_location, 
	ct.visualize = ct.visualize,                 ### selected cell types to visualize
	colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
	NumCols = 4,                                 ### number of columns in the figure panel
    pointSize = 3.0)                         ### point size in ggplot2 scatterplot  


#  visualize the cell type proportion correlation
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL)
p4
pdf("correlation.pdf")
p4
dev.off()

# imputation on the newly grided spatial locations
CARD_obj = CARD.imputation(CARD_obj,NumGrids = 6000,ineibor = 10,exclude = NULL)

location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
	y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
p5 <- ggplot(location_imputation, 
       aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    legend.position="bottom",
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))

# visualize the cell type proportion at an enhanced resolution
p6 <- CARD.visualize.prop(
	proportion = CARD_obj@refined_prop,                         
	spatial_location = location_imputation,            
	ct.visualize = ct.visualize,                    
	colors = c("lightblue","lightyellow","red"),    
	NumCols = 4)  

# visualize the marker gene expression at an enhanced resolution
p7 <- CARD.visualize.gene(
	spatial_expression = CARD_obj@refined_expression,
	spatial_location = location_imputation,
	gene.visualize = c("RELN","RXFP2","PDGFRA","KET5","KDR","ITGA7","MKI67","MRC1"),
	colors = NULL,
	NumCols = 6)
