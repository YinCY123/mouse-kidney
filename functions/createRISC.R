createRISC <- function(sce_list, ...){
    require(RISC)
    
    for(i in seq_along(sce_list)){
        current <- readsc(count = counts(sce_list[[i]]), 
                          cell = colData(sce_list[[i]]), 
                          gene = rowData(sce_list[[i]]), 
                          is.filter = F)
        # current@assay$logcounts <- logcounts(sce_list[[i]])
        sce_list[[i]] <- current
    }
    return(sce_list)
}
