SCE2cds <- function(SCE, ...){
    require(monocle3)
    require(SingleCellExperiment)
    require(dplyr)
    cds <- new_cell_data_set(expression_data = counts(SCE), 
                             cell_metadata = colData(SCE), 
                             gene_metadata = rowData(SCE) %>% as.data.frame %>% rename("gene_short_name" = symbol))
    
    # normalizaed data
    sn <- assayNames(SCE)
    for(name in sn[-1]){
        assay(cds, name) <- assay(SCE, name)
    }
    
    # reduced Dimensional data
    reduced_dims <- reducedDimNames(SCE)
    for(name in reduced_dims){
        reducedDim(cds, name) <- reducedDim(SCE, name)
    }
    return(cds)
}
