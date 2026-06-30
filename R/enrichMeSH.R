##' MeSH term enrichment analysis
##'
##'
##' @title enrichMeSH
##' @param gene a vector of entrez gene id
##' @param MeSHDb MeSHDb
##' @param database one of 'gendoo', 'gene2pubmed' or 'RBBH'
##' @param category one of "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L","M", "N", "V", "Z"
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param qvalueCutoff qvalue cutoff
##' @param meshdbVersion version of MeSH.db. If NULL(the default), use the latest version.
##' @return An \code{enrichResult} instance.
##' @importClassesFrom enrichit enrichResult
##' @importFrom enrichit ora_gson
##' @export
##' @seealso \code{class?enrichResult}
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
##' de <- names(geneList)[1:100]
##' x <- enrichMeSH(de, MeSHDb = db, database='gendoo', category = 'C')
##' }
##' @author Guangchuang Yu
enrichMeSH <- function(gene,
                       MeSHDb,
                       database = 'gendoo',
                       category = 'C',
                       pvalueCutoff=0.05,
                       pAdjustMethod="BH",
                       universe = NULL,
                       minGSSize = 10,
                       maxGSSize = 500,
                       qvalueCutoff = 0.2,
                       meshdbVersion = NULL) {

    MeSH_DATA <- get_MeSH_data(MeSHDb, database, category)

    res <- ora_gson(gene = gene,
                    pvalueCutoff=pvalueCutoff,
                    pAdjustMethod=pAdjustMethod,
                    universe = universe,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    qvalueCutoff = qvalueCutoff,
                    gson = MeSH_DATA
                    )

    return(res)
}


##' @importFrom yulab.utils get_fun_from_pkg
##' @importFrom gson gson
get_MeSH_data <- function(MeSHDb, database, category) {
    .meshesenv <- get_mesh_env()
    check_MeSHDb(MeSHDb)
    
    if (exists("meshtable", envir=.meshesenv)) {
        mesh <- get("meshtable", envir = .meshesenv)
    } else {
        ## db <- get_fun_from_pkg(MeSHDb, MeSHDb)
        mesh <- select(MeSHDb, keys=database, columns = c("GENEID", "MESHID","MESHCATEGORY"), keytype = "SOURCEDB")
        assign("meshtable", mesh, envir = .meshesenv)
    }

    category <- toupper(category)
    categories <- c("A", "B", "C", "D",
                    "E", "F", "G", "H",
                    "I", "J", "K", "L",
                    "M", "N", "V", "Z")

    if (!all(category %in% categories)) {
        stop("please check your 'category' parameter...")
    }

    mesh <- mesh[ mesh[,3] %in% category, ]
    mesh2gene <- mesh[, c(2,1)]

    ## meshdb <- get_fun_from_pkg("MeSH.db", "MeSH.db")
    mesh2name <- select(MeSHDb, keys=unique(mesh2gene[,1]), columns=c('MESHID', 'MESHTERM'), keytype='MESHID')

    gson(gsid2gene = mesh2gene, 
        gsid2name = mesh2name, 
        species = get_organism(MeSHDb),
        gsname = "MeSH")
}
