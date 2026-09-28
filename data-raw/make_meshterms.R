## Regenerate data/meshterms.rda from the annual NLM MeSH descriptor release.
## Update the URL each year, then review the row count and package checks.

url <- "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.gz"
tmp <- tempfile(fileext = ".gz")
download.file(url, tmp, mode = "wb", quiet = TRUE)
xml <- xml2::read_xml(gzfile(tmp))
records <- xml2::xml_find_all(xml, ".//DescriptorRecord")
meshterms <- data.frame(
    MESHID = xml2::xml_text(xml2::xml_find_first(records, "./DescriptorUI")),
    MESHTERM = xml2::xml_text(xml2::xml_find_first(records, "./DescriptorName/String")),
    stringsAsFactors = FALSE
)
meshterms <- meshterms[!is.na(meshterms$MESHID) & !is.na(meshterms$MESHTERM), ]
stopifnot(nrow(meshterms) > 30000, !anyDuplicated(meshterms$MESHID))
save(meshterms, file = "data/meshterms.rda", compress = "xz")
