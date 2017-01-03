##' semantic similarity between two gene vector
##'
##'
##' @title geneSim
##' @param geneID1 gene ID vector
##' @param geneID2 gene ID vector
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang" and "Lin"
##' @param combine One of "max", "avg", "rcmax", "BMA" methods, for combining semantic similarity scores of multiple DO terms associated with gene/protein.
##' @param semData gene annotation data for semantic measurement
##' @importFrom GOSemSim combineScores
##' @return score matrix
##' @export
##' @examples
##' ## hsamd <- meshdata("MeSH.Hsa.eg.db", category='A', computeIC=T, database="gendoo")
##' data(hsamd)
##' geneSim("241", "251", semData=hsamd, measure="Wang", combine="BMA")
##' @author Guangchuang Yu
geneSim <- function(geneID1,
                    geneID2=NULL,
                    measure="Wang",
                    combine="BMA",
                    semData) {

    if (!exists(".meshesEnv")) .initial()

    meshid1 <- lapply(geneID1, gene2MeSH, semData=semData)
    if (is.null(geneID2)) {
        geneID2 <- geneID1
        meshid2 <- meshid1
    } else {
        meshid2 <- lapply(geneID2, gene2MeSH, semData=semData)
    }

    m <- length(geneID1)
    n <- length(geneID2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- geneID1
    colnames(scores) <- geneID2

    for (i in 1:m) {
        if (length(geneID1) == length(geneID2) && all(geneID1 == geneID2)) {
           nn <- i
           flag <- TRUE
        } else {
            flag <- FALSE
            nn <- n
        }

        for (j in 1:nn) {
            if(any(!is.na(meshid1[[i]])) &&  any(!is.na(meshid2[[j]]))) {
                s <- meshSim(meshid1[[i]],
                           meshid2[[j]],
                           measure,
                           semData
                           )
                scores[i,j] = combineScores(s, combine)
                if (flag == TRUE && j != i) {
                    scores[j, i] <- scores[i,j]
                }
            }
        }
    }
    if (nrow(scores) == 1 & ncol(scores) == 1)
        scores = as.numeric(scores)
    return(scores)
}


gene2MeSH <- function(geneID, semData) {
    meshAnno <- semData@geneAnno
    meshAnno[meshAnno$GENEID == geneID, "MESHID"]
}
