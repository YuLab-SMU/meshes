##' semantic similarity between two MeSH term vectors
##'
##' 
##' @title meshSim
##' @param meshID1 MeSH term vector 
##' @param meshID2 MeSH term vector 
##' @param measure one of "Wang", "Resnik", "Rel", "Jiang" and "Lin"
##' @param semData annotation data for semantic measurement, output by meshdata function 
##' @return score matrix
##' @importFrom GOSemSim termSim
##' @export
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
meshSim <- function(meshID1,
                    meshID2,
                    measure = "Wang",
                    semData) {
    scores <- termSim(meshID1, meshID2, semData, measure)
    if(nrow(scores) == 1 & ncol(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}




## currently work
##
## library(GOSemSim)
## library(meshsim)
## hsamd = meshdata("MeSH.Hsa.eg.db", category='A', computeIC=T, database="gendoo")
## save(hsamd, file='hsamd.rda', compress='xz')
## meshSim("D000009", "D009130", semData=hsamd, measure="Resnik")
