setwd("/home/yincy/git/mouse-kidney/")

library(Seurat)
library(magrittr)
library(ggplot2)
kidney_integrated <- readRDS("../data/kidney/kidney-map/integrated_data/Seurat_integrated.rds")


for(i in c("cell_type", "location", "mito_percent")){
  DimPlot(kidney_integrated, 
          reduction = "tsne", 
          label = ifelse(i == "mito_percent", FALSE, TRUE), 
          group.by = i,
          label.size = 2) +
    NoLegend()
  
  ggsave(paste("figures/Seurat_", i, ".pdf", sep = ""))
}

