meshterms <- NULL
utils::data("meshterms", package = "meshes", envir = environment())
stopifnot(
    is.data.frame(meshterms),
    identical(names(meshterms), c("MESHID", "MESHTERM")),
    nrow(meshterms) > 30000,
    !anyDuplicated(meshterms$MESHID),
    all(nzchar(meshterms$MESHID)),
    all(nzchar(meshterms$MESHTERM))
)
