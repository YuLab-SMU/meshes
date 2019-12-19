##' Gene Set Enrichment Analysis of MeSH
##'
##'
##' @title gseMeSH
##' @param geneList order ranked geneList
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @param exponent weight of each step
##' @param nPerm permutation numbers
##' @param minGSSize minimal size of each geneSet for analyzing
##' @param maxGSSize maximal size of genes annotated for testing
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed logical
##' @param by one of 'fgsea' or 'DOSE'
##' @importClassesFrom DOSE gseaResult
##' @export
##' @return gseaResult object
##' @examples
##' \dontrun{
##' data(geneList, package="DOSE")
##' y <- gseMeSH(geneList, MeSHDb = "MeSH.Hsa.eg.db", database = 'gene2pubmed', category = "G")
##' }
##' @author Yu Guangchuang
gseMeSH <- function(geneList,
                    MeSHDb,
                    database = 'gendoo',
                    category = 'C',
                    exponent      = 1,
                    nPerm         = 1000,
                    minGSSize     = 10,
                    maxGSSize     = 500,
                    pvalueCutoff=0.05,
                    pAdjustMethod="BH",
                    verbose       = TRUE,
                    seed          = FALSE,
                    by = 'fgsea') {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)
    res <-  GSEA_internal(geneList = geneList,
                          exponent = exponent,
                          nPerm = nPerm,
                          minGSSize = minGSSize,
                          maxGSSize = maxGSSize,
                          pvalueCutoff = pvalueCutoff,
                          pAdjustMethod = pAdjustMethod,
                          verbose = verbose,
                          USER_DATA = MeSH_DATA,
                          seed = seed,
                          by = by)

    meshdb <- get_fun_from_pkg("MeSH.db", "MeSH.db")
    id <- res@result$ID
    mesh2name <- select(meshdb, keys=id, columns=c('MESHID', 'MESHTERM'), keytype='MESHID')
    res@result$Description <- mesh2name[match(id, mesh2name[,1]), 2]
    res@organism <- get_organism(MeSHDb)
    res@setType <- "MeSH"

    return(res)
}
