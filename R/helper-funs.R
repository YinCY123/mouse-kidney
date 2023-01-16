# convert between symbol id and ensembl id

idconv <- function(from_id, from_type = "SYMBOL", to_type = "GENEID", db = EnsDb.Mmusculus.v79, ...){
    require(EnsDb.Mmusculus.v79)
    require(magrittr)
    
    mapIds(x = db, 
           keys = from_id, 
           keytype = from_type, 
           column = to_type) %>% 
        unname()
}
