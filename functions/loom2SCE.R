list2df <- function(list){ 
    maxL <- sapply(list, length) %>% max()
    for(i in seq_along(list)){
        l = length(list[[i]])
        if(l < maxL){
            list[[i]] <- append(list[[i]], values = rep(NA, (maxL - l)))
        }
    }
    list <- as.data.frame(list)
    # colnames(list) <- NULL
    # list <- t(list)
    # colnames(list) <- paste(rep("order", times = maxL), 0:(maxL - 1), sep = "_")
    return(as.data.frame(list))
}


loom2SCE <- function(loom, ...){
    require(SingleCellExperiment)
    require(loomR)
    require(magrittr)
    require(Matrix)
    
    col_attrs_names_ids <- loom[["col_attrs"]] %>% as.list() %>% .[["names"]] %>% grep("^[^_]", .)
    col_attrs_names <- loom[["col_attrs"]] %>% as.list() %>% .[["names"]] %>% .[col_attrs_names_ids]
    
    row_attrs_names_ids <- loom[["row_attrs"]] %>% as.list() %>% .[["names"]] %>% grep("^[^_]", .)
    row_attrs_names <- loom[["row_attrs"]] %>% as.list() %>% .[["names"]] %>% .[row_attrs_names_ids]
    
    col_attrs_df <- list()
    for(i in 1:length(col_attrs_names)){
        col_attrs_df[[col_attrs_names[[i]]]] <- loom[[paste("col_attrs", col_attrs_names[[i]], sep = "/")]][]
    }
    
    if(sum(abs(diff(sapply(col_attrs_df, length)))) != 0){
        col_attrs_df <- list2df(col_attrs_df)
    }else{
        col_attrs_df <- as.data.frame(col_attrs_df)
    }
    rownames(col_attrs_df) <- paste(col_attrs_df$CellID, loom[["col_attrs/emptydrops_FDR"]][], sep = "_")
    print("done...")
    
    row_attrs_df <- list()
    for(i in 1:length(row_attrs_names)){
        row_attrs_df[[row_attrs_names[[i]]]] <- loom[[paste("row_attrs", row_attrs_names[[i]], sep = "/")]][]
    }
    
    if(sum(abs(diff(sapply(row_attrs_df, length)))) != 0){
        row_attrs_df <- list2df(row_attrs_df)
    }else{
        row_attrs_df <- as.data.frame(row_attrs_df)
    }
    rownames(row_attrs_df) <- row_attrs_df$Accession
    
    
    counts <- t(loom[["matrix"]][,])
    colnames(counts) <- paste(loom[["col_attrs/CellID"]][], loom[["col_attrs/emptydrops_FDR"]][], sep = "_")
    rownames(counts) <- loom[["row_attrs/Accession"]][]
    
    sce <- SingleCellExperiment(assays = list(counts = as(counts, "CsparseMatrix")), 
                                rowData = row_attrs_df, 
                                colData = col_attrs_df)
    return(sce)
}




