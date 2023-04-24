harmonizeSCE <- function(sce_list, method, ...){
    require(SingleCellExperiment)
    require(Matrix)
    require(AnnotationHub)
    ensdb <- AnnotationHub()[["AH104895"]]
    
    ugenes <- Reduce(union, lapply(sce_list, rownames))
    
    # for harmony
    if(method == "harmony"){
        counts <- list()
        # equal genes
        for(i in seq_along(sce_list)){
            dif_genes <- setdiff(ugenes, rownames(sce_list[[i]]))
            dif_mat <- matrix(data = 0, 
                              nrow = length(dif_genes), 
                              ncol = length(colnames(sce_list[[i]])), 
                              dimnames = list(dif_genes, colnames(sce_list[[i]])))
            dif_mat <- as(dif_mat, "CsparseMatrix")
            
            # assay
            counts[[i]] <- rbind(counts(sce_list[[i]]), dif_mat)
        }
        return(Reduce(cbind, counts))  
    }
    # return(Reduce(cbind, counts))
    
    # for RISC
    if(method == "RISC"){
        for(i in seq_along(sce_list)){
            dif_gene <- setdiff(ugenes, rownames(sce_list[[i]]))
            
            dif_mat <- matrix(data = 0, 
                              nrow = length(dif_gene), 
                              ncol = length(colnames(sce_list[[i]])), 
                              dimnames = list(dif_gene, colnames(sce_list[[i]])))
            dif_mat <- as(dif_mat, "CsparseMatrix")
            counts <- rbind(counts(sce_list[[i]]), dif_mat)
            # logcounts <- rbind(logcounts(sce_list[[i]]), dif_mat)
            
            gene <- data.frame(ensembl = rownames(counts), 
                               symbol = mapIds(ensdb, keys = rownames(counts), keytype = "GENEID", column = "SYMBOL"),
                               row.names = rownames(counts))
            
            cell <- data.frame(cell_ids = colnames(counts),
                               cell_type = colData(sce_list[[i]])$cell_type,
                               row.names = colnames(counts))
            
            sce_list[[i]] <- SingleCellExperiment(assays = list(counts = counts), 
                                                  rowData = gene, 
                                                  colData = cell)
        }
        
        # reorder
        for(i in seq_along(sce_list)){
            ids <- match(rownames(sce_list[[1]]), rownames(sce_list[[i]]))
            sce_list[[i]] <- sce_list[[i]][ids, ]
        }
        return(sce_list)
    }
    
    # for fastMNN
    # TODO
}
