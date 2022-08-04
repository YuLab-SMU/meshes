##' download the latest version of mesh data and stored in a 'GSON' object
##'
##'
##' @title gson_Reactome
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @importFrom utils stack
##' @importFrom gson gson
##' @export
gson_mesh <- function(MeSHDb, database, category) {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category) 
    EXTID2PATHID = get("EXTID2PATHID", envir = MeSH_DATA)
    PATHID2EXTID = get("PATHID2EXTID", envir = MeSH_DATA)
    PATHID2NAME = get("PATHID2NAME", envir = MeSH_DATA) # NULL

    reactomeAnno <- stack(PATHID2EXTID)   
    gsid2gene <- reactomeAnno[, c(2,1)]
    colnames(gsid2gene) <- c("gsid", "gene")
    gsid2gene <- unique(gsid2gene[!is.na(gsid2gene[,2]), ]) 

    termname <- PATHID2NAME
    gsid2name <- NULL
    if (is.null(termname)) {
        gsid2name <- data.frame(gsid = names(termname), name = termname)
    }

    gson(gsid2gene = gsid2gene, 
        gsid2name = gsid2name,
        # species = species,
        gsname = "MeSH",
        # version = version,
        accessed_date = as.character(Sys.Date())
        # keytype = keytype
    )
}