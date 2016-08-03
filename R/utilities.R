.initial <- function() {
    assign(".meshsimEnv", new.env(), .GlobalEnv)
    tryCatch(utils::data(list="meshtbl",
                         package = "meshsim"))
    meshtbl <- get("meshtbl")
    assign("meshtbl", meshtbl, envir = .meshsimEnv)
    rm(meshtbl, envir=.GlobalEnv)
}

