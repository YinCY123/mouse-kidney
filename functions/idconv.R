# convert between symbol id and ensembl id

idconv <- function(from_id, from_type = "SYMBOL", to_type = "GENEID", db = EnsDb.Mmusculus.v79, ...){
    require(AnnotationHub)
    require(magrittr)
    ensdb <- AnnotationHub(localHub = TRUE)[["AH104895"]]
    
    mapIds(x = ensdb, 
           keys = from_id, 
           keytype = from_type, 
           column = to_type) %>% 
        unname()
}
