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
##' @export
##' @examples
##' \dontrun{
##' library(meshes)
##' library(AnnotationHub)
##' ah <- AnnotationHub()
##' qr_hsa <- query(ah, c("MeSHDb", "Homo sapiens"))
##' ## inspect qr_hsa and select the organism-specific MeSHDb resource
##' filepath_hsa <- qr_hsa[[1]]
##' db <- MeSHDbi::MeSHDb(filepath_hsa)
##' hsamd <- meshdata(db, category='A', computeIC=T, database="gendoo")
##' }
##' @author Guangchuang Yu 
meshdata <- function(MeSHDb=NULL, database, category, computeIC = FALSE) {
    if (is.null(MeSHDb)) {
        return(new("GOSemSimDATA",
                   ont = category))
    }

    check_MeSHDb(MeSHDb)

    # MeSHDb <- load_OrgDb(MeSHDb)
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

check_MeSHDb <- function(MeSHDb) {
    expected <- c("GENEID", "MESHCATEGORY", "MESHID", "SOURCEDB")
    keytypes <- tryCatch(
        AnnotationDbi::keytypes(MeSHDb),
        error = function(e) e
    )

    if (inherits(keytypes, "error") || !all(expected %in% keytypes)) {
        msg <- paste0(
            "`MeSHDb` should be an organism-specific MeSHDb annotation ",
            "database containing a DATA table with columns: ",
            paste(expected, collapse = ", "), ". ",
            "When using AnnotationHub, inspect the query result and select ",
            "the organism-specific MeSHDb resource instead of assuming the ",
            "first result is suitable."
        )
        if (inherits(keytypes, "error")) {
            msg <- paste0(msg, " Original error: ", conditionMessage(keytypes))
        }
        stop(msg, call. = FALSE)
    }

    invisible(TRUE)
}

getOffsprings <- function(meshID) {
    children <- get_mesh_children()
    res <- c()
    id <- meshID
    while(any(id %in% names(children))) {
        cid <- unlist(children[id], use.names = FALSE)
        res <- c(res, cid)
        id <- cid
    }
    return(unique(res))
}

getAncestors <- function(meshID) {
    parents <- get_mesh_parents()
    res <- c()
    id <- meshID
    while(any(id %in% names(parents))) {
        pid <- unlist(parents[id], use.names = FALSE)
        res <- c(res, pid)
        id <- pid
    }
    return(unique(res))
}

get_mesh_children <- function() {
    .meshesEnv <- get_mesh_env()
    if (!exists("mesh_children", envir = .meshesEnv, inherits = FALSE)) {
        meshtbl <- get("meshtbl", envir = .meshesEnv)
        assign("mesh_children", split(meshtbl$meshID, meshtbl$parent), envir = .meshesEnv)
    }
    get("mesh_children", envir = .meshesEnv)
}

get_mesh_parents <- function() {
    .meshesEnv <- get_mesh_env()
    if (!exists("mesh_parents", envir = .meshesEnv, inherits = FALSE)) {
        meshtbl <- get("meshtbl", envir = .meshesEnv)
        assign("mesh_parents", split(meshtbl$parent, meshtbl$meshID), envir = .meshesEnv)
    }
    get("mesh_parents", envir = .meshesEnv)
}



get_mesh_term_names <- function(MeSHDb, mesh_ids) {
    mesh_ids <- unique(as.character(mesh_ids))
    if (!length(mesh_ids)) {
        return(data.frame(MESHID = character(), MESHTERM = character(),
                          stringsAsFactors = FALSE))
    }

    available <- tryCatch(
        AnnotationDbi::columns(MeSHDb),
        error = function(e) character()
    )

    if ("MESHTERM" %in% available) {
        res <- AnnotationDbi::select(
            MeSHDb,
            keys = mesh_ids,
            columns = c("MESHID", "MESHTERM"),
            keytype = "MESHID"
        )
        res <- res[!duplicated(res$MESHID), c("MESHID", "MESHTERM"),
                   drop = FALSE]
        res$MESHID <- as.character(res$MESHID)
        res$MESHTERM <- as.character(res$MESHTERM)
        return(res[match(mesh_ids, res$MESHID), , drop = FALSE])
    }

    ## Since the current AnnotationHub MeSHDb schema contains only the
    ## gene-to-MeSH correspondence, use the versioned NLM descriptor map
    ## bundled with meshes for stable, offline term names.
    .meshesEnv <- get_mesh_env()
    if (!exists("meshterms", envir = .meshesEnv, inherits = FALSE)) {
        utils::data("meshterms", package = "meshes", envir = .meshesEnv)
    }
    meshterms <- get("meshterms", envir = .meshesEnv, inherits = FALSE)
    res <- meshterms[match(mesh_ids, meshterms$MESHID),
                     c("MESHID", "MESHTERM"), drop = FALSE]
    res$MESHID <- mesh_ids
    missing <- is.na(res$MESHTERM) | !nzchar(res$MESHTERM)
    if (any(missing)) {
        ## Preserve a usable gene-set object for future descriptor additions
        ## that are newer than the bundled NLM release.
        res$MESHTERM[missing] <- res$MESHID[missing]
        warning(sum(missing), " MeSH IDs have no bundled descriptor name; ",
                "using the MeSH ID as the term label.", call. = FALSE)
    }
    res
}


#' Get the current organism-specific MeSHDb from AnnotationHub
#'
#' @param meshdbVersion optional AnnotationHub release/tag to select
#' @noRd
get_meshdb <- function(meshdbVersion = NULL) {
    .meshesEnv <- get_mesh_env()
    if (exists("meshdb", envir = .meshesEnv, inherits = FALSE) &&
        exists("meshdbVersion", envir = .meshesEnv, inherits = FALSE)) {
        meshdbVersion2 <- get("meshdbVersion", envir = .meshesEnv)
        if (identical(meshdbVersion, meshdbVersion2)) {
            return(get("meshdb", envir = .meshesEnv))
        }
    }

    ah <- AnnotationHub::AnnotationHub()
    query_terms <- c("MeSHDb", "Homo sapiens")
    if (!is.null(meshdbVersion)) {
        query_terms <- c(query_terms, meshdbVersion)
    }
    records <- AnnotationHub::query(ah, query_terms)
    if (!length(records)) {
        stop("No human MeSHDb record matched the requested release.",
             call. = FALSE)
    }

    dates <- as.Date(records$rdatadateadded)
    record_index <- if (all(is.na(dates))) length(records) else which.max(dates)
    meshdb <- MeSHDbi::MeSHDb(records[[record_index]])
    assign("meshdb", meshdb, envir = .meshesEnv)
    assign("meshdbVersion", meshdbVersion, envir = .meshesEnv)
    meshdb
}


get_mesh_env <- function () {
    if (!exists(".meshesEnv")) .initial()
    get(".meshesEnv", envir = .GlobalEnv)
}


