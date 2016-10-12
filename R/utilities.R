.initial <- function() {
    pos <- 1
    envir <- as.environment(pos) 
    
    assign(".meshesEnv", new.env(), envir = envir)
    .meshesEnv <- get(".meshesEnv", envir = .GlobalEnv)
    tryCatch(utils::data(list="meshtbl",
                         package = "meshes"))
    meshtbl <- get("meshtbl")
    assign("meshtbl", meshtbl, envir = .meshesEnv)
    rm(meshtbl, envir=.GlobalEnv)
}

##' construct annoData for semantic measurement
##'
##' 
##' @title meshdata
##' @param MeSHDb MeSHDb package 
##' @param database one of supported database
##' @param category one of supported category
##' @param computeIC logical value
##' @return a GOSemSimDATA object
##' @importFrom AnnotationDbi metadata
##' @importFrom AnnotationDbi keys
##' @importFrom methods new
##' @importClassesFrom GOSemSim GOSemSimDATA
##' @importFrom GOSemSim load_OrgDb
##' @export
##' @examples
##' meshdata("MeSH.Cel.eg.db", category='A', computeIC=FALSE, database="gene2pubmed")
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

    meshAnno <- select(MeSHDb, keys=kk, keytype="GENEID", columns = c("GENEID", "MESHCATEGORY", "MESHID", "SOURCEDB"))
    meshAnno <- meshAnno[meshAnno$SOURCEDB %in% database,]
    meshAnno <- meshAnno[meshAnno$MESHCATEGORY == category, ]

    res <- new("GOSemSimDATA",
               keys = kk,
               ont = category,
               geneAnno = meshAnno,
               metadata = metadata(MeSHDb)
               )
    
    if (computeIC) {
        res@IC = computeIC(meshAnno, category)
    }
    
    return(res)
}

## @importFrom MeSH.AOR.db MeSH.AOR.db
##' @importFrom AnnotationDbi select
computeIC <- function(meshAnno, category) {
    meshdata <- get("meshtbl", envir=.meshesEnv)
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
    meshtbl <- get("meshtbl", envir=.meshesEnv)
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
    meshtbl <- get("meshtbl", envir=.meshesEnv)
    res <- c()
    id <- meshID
    while(any(id %in% meshtbl$meshID)) {
        pid <- meshtbl[meshtbl$meshID %in% id, "parent"]
        res <- c(res, pid)
        id <- pid
    }
    return(unique(res))
}


