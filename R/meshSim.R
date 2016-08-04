meshSim <- function(meshID1,
                    meshID2,
                    measure = "Wang",
                    meshdata) {
    scores <- termSim(meshID1, meshID2, meshdata, measure)
    if(nrow(scores) == 1 & ncol(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}

termSim <- GOSemSim::termSim


## currently work
##
## library(GOSemSim)
## library(meshsim)
## md = meshdata("MeSH.Hsa.eg.db", category='A', computeIC=T, database="gene2pubmed")
## meshSim("D000009", "D009130", meshdata=md, measure="Resnik")
