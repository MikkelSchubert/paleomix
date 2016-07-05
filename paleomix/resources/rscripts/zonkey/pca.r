args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    cat("Usage: plot_pca.R <input_prefix> <name_table> <output_file>\n",
        file=stderr())
    quit(status=1)
}


library(ggplot2)
library(ggrepel)


expanded.range <- function(values, expand.by=0.2)
{
    r <- range(values)
    p <- abs(r[2] - r[1]) * expand.by

    return(r + c(-0.5, 0.5) * p)
}


calc.lims <- function(x, y)
{
    x <- expanded.range(x)
    y <- expanded.range(y)
    xlen <- x[2] - x[1]
    ylen <- y[2] - y[1]

    if (xlen < ylen) {
        x <- x + c(-0.5, 0.5) * (ylen - xlen)
    }

    return(x)
}


plot.pca <- function(input_prefix, names.table=NULL)
{
    d <- read.table(paste(input_prefix, "evec", sep="."), as.is=TRUE)
    v <- read.table(paste(input_prefix, "eval", sep="."), as.is=TRUE)

    pc1 <- round(v$V1[1] / sum(v$V1), 3)
    pc2 <- round(v$V1[2] / sum(v$V1), 3)

    d$V1 <- sapply(strsplit(d$V1, ":", fixed=TRUE), function(x) x[1])
    d <- d[, 1:3]
    colnames(d) <- c("Name", "X", "Y")

    colors <- read.table(names.table, comment.char="", header=TRUE)
    final <- merge(d, colors)

    pp <- ggplot(final)
    pp <- pp + geom_text_repel(aes(x=X, y=Y, label=Name))

    pp <- pp + geom_point(aes(x=X, y=Y), color=final$Color, size=3)

    pp <- pp + xlab(sprintf("PC1: %.1f%%", pc1 * 100))
    pp <- pp + ylab(sprintf("PC2: %.1f%%", pc2 * 100))

    pp <- pp + theme_minimal()
    pp <- pp + xlim(calc.lims(final$X, final$Y))
    pp <- pp + ylim(calc.lims(final$Y, final$X))

    # Force 1:1 aspect ratio for the PCA plot; this must be done in
    # order to prevent misrepresenting distances on different axes.
    pp <- pp + coord_fixed()

    return(pp)
}


input_prefix <- args[1]
names_table <- args[2]
output_prefix <- args[3]

pdf(file=paste(output_prefix, ".pdf", sep=""))
plot.pca(input_prefix, names_table)
dev.off()

# bitmap is preferred, since it works in a headless environment
bitmap(paste(output_prefix, ".png", sep=""), height=6, width=6, res=96, taa=4, gaa=4)
plot.pca(input_prefix, names_table)
dev.off()
