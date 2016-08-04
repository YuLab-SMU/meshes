.initial <- function() {
    assign(".meshsimEnv", new.env(), .GlobalEnv)
    tryCatch(utils::data(list="meshtbl",
                         package = "meshsim"))
    meshtbl <- get("meshtbl")
    assign("meshtbl", meshtbl, envir = .meshsimEnv)
    rm(meshtbl, envir=.GlobalEnv)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title meshdata
##' @param MeSHDb MeSHDb package 
##' @param database one of supported database
##' @param category one of supported category
##' @param computeIC logical value
##' @return a GOSemSimDATA object
##' @importFrom AnnotationDbi metadata
##' @importFrom AnnotationDbi keys
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu 
meshdata <- function(MeSHDb=NULL, database, category, computeIC = FALSE) {
    if (is.null(MeSHDb)) {
        return(new("GOSemSimDATA",
                   ont = category))
    }

    MeSHDb <- load_OrgDb(MeSHDb)
    SOURCEDB <- keys(MeSHDb, keytype="SOURCEDB")
    if (!database %in% SOURCEDB) {
        msg <- paste0("supported database is/are '", paste(SOURCEDB, sep='/'), "', input parameter not matched...")
        stop(msg)
    }
    
    kk <- as.character(keys(MeSHDb, keytype="GENEID"))

    meshAnno <- select(MeSHDb, keys=kk, keytype="GENEID", columns = c("MESHCATEGORY", "MESHID", "SOURCEDB"))
    meshAnno <- meshAnno[meshAnno$SOURCEDB %in% database,]
    meshAnno <- meshAnno[meshAnno$MESHCATEGORY == category, ]

    res <- new("GOSemSimDATA",
               keys = kk,
               ont = category,
               metadata = metadata(MeSHDb))
    if (computeIC) {
        res@IC = computeIC(meshAnno, category)
    }

    return(res)
}

##' @importFrom MeSH.AOR.db MeSH.AOR.db
##' @importFrom AnnotationDbi select
computeIC <- function(meshAnno, category) {
    meshdata <- get("meshtbl", envir=.meshsimEnv)
    meshids <- unique(meshdata[meshdata$Ontology == category, "meshID"])

    meshterms <- meshAnno$MESHID
    meshcount <- table(meshterms)
    meshname <- names(meshcount)

    mesh.diff <- setdiff(meshids, meshname)
    m <- double(length(mesh.diff))
    names(m) <- mesh.diff
    meshcount <- as.vector(meshcount)
    names(meshcount) <- meshname

    ## offspring.df <- select(MeSH.AOR.db, keys=category, columns=c("ANCESTOR", "OFFSPRING"), keytype="CATEGORY")
    ## Offsprings <- split(offspring.df$OFFSPRING, offspring.df$ANCESTOR) 
    Offsprings <- lapply(meshids, getOffsprings)
    names(Offsprings) <- meshids
    cnt <- meshcount[meshids] + sapply(meshids, function(i) sum(meshcount[Offsprings[[i]]], na.rm=TRUE))
    names(cnt) <- meshids

    p <- cnt/sum(meshcount)
    IC <- -log(p)
    return(IC)
}

getOffsprings <- function(meshID) {
    meshtbl <- get("meshtbl", envir=.meshsimEnv)
    res <- c()
    id <- meshID
    while(any(id %in% meshtbl$parent)) {
        cid <- meshtbl[meshtbl$parent %in% id, "meshID"]
        res <- c(res, cid)
        id <- cid
    }
    return(unique(res))
}

getAncestors <- function(meshID) {
    meshtbl <- get("meshtbl", envir=.meshsimEnv)
    res <- c()
    id <- meshID
    while(any(id %in% meshtbl$meshID)) {
        pid <- meshtbl[meshtbl$meshID %in% id, "parent"]
        res <- c(res, pid)
        id <- pid
    }
    return(unique(res))
}

load_OrgDb <- GOSemSim:::load_OrgDb
