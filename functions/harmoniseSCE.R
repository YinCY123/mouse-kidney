harmonizeSCE <- function(sce_list, method, ...){
    require(SingleCellExperiment)
    require(Matrix)
    require(EnsDb.Mmusculus.v79)
    # require(AnnotationHub)
    # ensdb <- AnnotationHub(localHub = FALSE)[["AH104895"]]
    
    ugenes <- Reduce(union, lapply(sce_list, rownames)) %>% unique()
    
    
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
            
            gene <- data.frame(ensembl = rownames(counts), 
                               symbol = mapIds(EnsDb.Mmusculus.v79, keys = rownames(counts), keytype = "GENEID", column = "SYMBOL"),
                               row.names = rownames(counts))
            
            cell <- data.frame(cell_ids = colnames(counts),
                               cell_type = colData(sce_list[[i]])$cell_type,
                               mito_percent = colData(sce_list[[i]])$mito_percent,
                               location = colData(sce_list[[i]])$location,
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
    if(method == "fastMNN"){
        for(name in names(sce_list)){
            dif_gene <- setdiff(ugenes, rownames(sce_list[[name]]))
            
            dif_mat <- matrix(data = 0, 
                              nrow = length(dif_gene), 
                              ncol = length(colnames(sce_list[[name]])), 
                              dimnames = list(dif_gene, colnames(sce_list[[name]])))
            dif_mat <- as(dif_mat, "CsparseMatrix")
            counts <- rbind(counts(sce_list[[name]]), dif_mat)
            
            gene <- data.frame(ensembl = rownames(counts), 
                               symbol = mapIds(EnsDb.Mmusculus.v79, keys = rownames(counts), keytype = "GENEID", column = "SYMBOL"),
                               row.names = rownames(counts))
            
            geo_accession <- c("dev_cell" = "GSE129798", 
                               "glomer" = "GSE146912", 
                               "drg" = "GSE175421",
                               "sym" = "GSE78845")
            cell <- data.frame(cell_ids = colnames(counts),
                               cell_type = colData(sce_list[[name]])$cell_type,
                               mito_percent = colData(sce_list[[name]])$mito_percent,
                               location = colData(sce_list[[name]])$location,
                               geo_accession = geo_accession[[name]],
                               row.names = colnames(counts))
            
            sce_list[[name]] <- SingleCellExperiment(assays = list(counts = counts), 
                                                     rowData = gene, 
                                                     colData = cell)
        }
        
        # reorder
        for(name in names(sce_list)){
            ids <- match(rownames(sce_list[[1]]), rownames(sce_list[[name]]))
            sce_list[[name]] <- sce_list[[name]][ids, ]
        }
        return(sce_list)
    }
}


