##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    pkgVersion <- packageDescription(pkgname, fields="Version")
    msg <- paste0(pkgname, " v", pkgVersion, "  ", "\n\n")
    citation <- paste0("If you use ", pkgname, " in published research, please cite the most appropriate paper(s):\n\n",
                       "Guangchuang Yu.",
                       "Using meshes for MeSH term enrichment and semantic analyses.",
                       "Bioinformatics 2018, 34(21):3766-3767, doi:10.1093/bioinformatics/bty410")
    
    packageStartupMessage(paste0(msg, citation))
}

.onLoad <- function(libname, pkgname) {
    .initial()
} 
