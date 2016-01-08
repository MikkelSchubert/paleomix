args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    cat("Usage: requires.R <module>\n")
    quit(status=1)
}

for (module in args) {
    if (!require(module, character.only=TRUE)) {
        quit(status=1)
    }

    # Magic string to detect successful loading
    cat(paste("d0fd3ea6:", packageVersion(module), " \n"))

    quit(status=0)
}
