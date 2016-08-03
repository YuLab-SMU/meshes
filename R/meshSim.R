meshSim <- function(meshID1,
                    meshID2,
                    measure = "Wang",
                    category) {
    scores <- termSim(meshID1, meshID2, meshdata(category), measure)
    if(nrow(scores) == 1 & ncol(scores) == 1)
        scores <- as.numeric(scores)
    return(scores)
}

termSim <- GOSemSim::termSim

meshdata <- function(category) {
    new("GOSemSimDATA",
        ont = category)
}

## currently work
##
## library(GOSemSim)
## library(meshsim)
## meshSim("D004312", "D009852", category='A')
