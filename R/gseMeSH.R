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
##' @param eps boundary for calculating the p value in multilevel mode
##' @param pvalueCutoff pvalue Cutoff
##' @param pAdjustMethod pvalue adjustment method
##' @param verbose print message or not
##' @param seed random seed for reproducibility, set to a number (or TRUE to use a
##'   fixed default seed) to make the result reproducible, or FALSE (default) to draw
##'   a random seed on each run, so results may vary between runs. The underlying
##'   permutation engine uses its own RNG seeded with this value; see
##'   \code{enrichit::gsea()} for details.
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
##' ## inspect qr_hsa and select the organism-specific MeSHDb resource
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
                    eps           = 1e-10,
                    pvalueCutoff  = 0.05,
                    pAdjustMethod = "BH",
                    verbose       = TRUE,
                    seed          = FALSE,
                    ...) {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)
    
    res <-  gsea_gson(geneList         = geneList,
                      gson             = MeSH_DATA,
                      nPerm            = nPerm,
                      exponent         = exponent,
                      minGSSize        = minGSSize,
                      maxGSSize        = maxGSSize,
                      eps              = eps,
                      pvalueCutoff     = pvalueCutoff,
                      pAdjustMethod    = pAdjustMethod,
                      verbose          = verbose,
                      seed             = seed,
                      ...)

    if (is.null(res))
        return(res)

    res@setType <- "MeSH"

    return(res)
}
