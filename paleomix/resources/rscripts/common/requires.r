args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    cat("Usage: requires.R <module>\n")
    quit(status=1)
}

packages <- installed.packages()

for (module in args) {
    if (!(module %in% packages)) {
        cat(paste("R module '", module, "' is not installed!\n", sep=''))
        cat(paste("Please run R and execute the following command:\n\n"))
        cat(paste("> install.packages('", module, "')\n", sep=''))
        quit(status=1)
    }

    # Magic string to detect successful loading
    cat(paste("d0fd3ea6:", packageVersion(module), " \n"))

    quit(status=0)
}
