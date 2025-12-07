##' Gene Set Enrichment Analysis of MeSH
##'
##'
##' @title gseMeSH
##' @param geneList order ranked geneList
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @param nPerm number of permutations.
##' @param exponent weight of each step
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param ... other parameter
##' @importClassesFrom enrichit gseaResult
##' @importFrom enrichit gsea_gson
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
##' }
##' @author Yu Guangchuang
gseMeSH <- function(geneList,
                    MeSHDb,
                    database      = 'gendoo',
                    category      = 'C',
                    nPerm         = 1000,
                    exponent      = 1,
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",
                    verbose       = TRUE,
                    ...) {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)
    
    res <-  gsea_gson(geneList         = geneList,
                      gson             = MeSH_DATA,
                      nPerm            = nPerm,
                      exponent         = exponent,
                      minGSSize        = minGSSize,
                      maxGSSize        = maxGSSize,
                      pvalueCutoff     = pvalueCutoff,
                      pAdjustMethod    = pAdjustMethod,
                      verbose          = verbose,
                      ...)
    
    res@setType <- "MeSH"

    return(res)
}
