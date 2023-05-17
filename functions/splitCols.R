splitCols <- function(sce, f, ...){
    v <- vector(mode = "list", length = length(unique(f)))
    f <- split(f, f)
    names(v) <- names(f)
    
    for(name in names(v)){
        v[[name]] <- sce[, as.numeric(f[[name]])]
    }
    
    return(v)
}