##' Gene Set Enrichment Analysis of MeSH
##'
##'
##' @title gseMeSH
##' @param geneList order ranked geneList
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param eps This parameter sets the boundary for calculating the p value.
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @param meshdbVersion version of MeSH.db. If NULL(the default), use the latest version.
##' @param gson a GSON object.
##' @param ... other parameter
##' @importClassesFrom DOSE gseaResult
##' @export
##' @return gseaResult object
##' @examples
##' \dontrun{
##' library(meshes)
##' library(AnnotationHub)
##' ah <- AnnotationHub()
##' qr_hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
##' filepath_hsa <- qr_hsa[[1]]
##' db <- MeSHDbi::MeSHDb(filepath_hsa)
##' data(geneList, package="DOSE")
##' y <- gseMeSH(geneList, MeSHDb = db, database = 'gene2pubmed', category = "G")

##' gsonMesh <- gson_mesh(MeSHDb = db, database = 'gene2pubmed', category = "G")
##' y2 <- gseMeSH(geneList, gson = gson)
##' }
##' @author Yu Guangchuang
gseMeSH <- function(geneList,
                    MeSHDb,
                    database      = 'gendoo',
                    category      = 'C',
                    exponent      = 1,
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    eps           = 1e-10,
                    pvalueCutoff  =0.05,
                    pAdjustMethod ="BH",
                    verbose       = TRUE,
                    seed          = FALSE,
                    by            = 'fgsea',
                    meshdbVersion = NULL,
                    gson = NULL,
                    ...) {
    if (is.null(gson)) {
        MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)      
        species <- get_organism(MeSHDb)  
    } else {
        MeSH_DATA <- gson
        species <- MeSH_DATA@sapiens
    }
   
    
    res <-  GSEA_internal(geneList         = geneList,
                          exponent         = exponent,
                          minGSSize        = minGSSize,
                          maxGSSize        = maxGSSize,
                          eps              = eps,
                          pvalueCutoff     = pvalueCutoff,
                          pAdjustMethod    = pAdjustMethod,
                          verbose          = verbose,
                          USER_DATA        = MeSH_DATA,
                          seed             = seed,
                          by               = by,
                          ...)
    

    # meshdb <- get_fun_from_pkg("MeSH.db", "MeSH.db")
    meshdb <- get_meshdb(meshdbVersion = meshdbVersion)
    id <- res@result$ID
    mesh2name <- select(meshdb, keys=id, columns=c('MESHID', 'MESHTERM'), keytype='MESHID')
    res@result$Description <- mesh2name[match(id, mesh2name[,1]), 2]
    res@organism <- species
    res@setType <- "MeSH"

    return(res)
}
