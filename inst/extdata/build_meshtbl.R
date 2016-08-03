library(MeSH.PCR.db)
kts <- keytypes(MeSH.PCR.db)
kt <- kts[2]

cls = columns(MeSH.PCR.db)

ks <- keys(MeSH.PCR.db, keytype=kts[2])
res <- select(MeSH.PCR.db, keys=ks, columns=cls, keytype=kt)

meshtbl <- meshtbl[, c(2, 3, 1)]
colnames(meshtbl) <- c("meshID", "parent", "Ontology")

save("meshtbl", file="meshtbl.rda", compress='xz')
