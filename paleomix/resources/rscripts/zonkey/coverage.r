#!/usr/bin/env Rscript
# Copyright (c) 2013 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    cat("Usage: test.R <table> <output_prefix>\n")
    quit(status=1)
}

library(ggplot2)


plot_coverage <- function(filename)
{
    tbl <- read.table(filename, as.is=TRUE, header=TRUE)
    tbl$Hits <- as.numeric(tbl$Hits)
    tbl$Size <- as.numeric(tbl$Size)
    tbl$Sample <- as.factor("Sample")

    # Correct size by number of uncalled ('N' / '-') bases
    tbl$RelHits <-  tbl$Hits / (tbl$Size - tbl$Ns)

    autosomes <- tbl[tbl$ID != 'X' & tbl$ID != 'Y',]
    autosomes$ID <- as.numeric(autosomes$ID)
    autosomes$NormHits <- 2 * autosomes$RelHits / mean(autosomes$RelHits)

    pp <- ggplot()

    pp <- pp + geom_hline(yintercept=range(autosomes$NormHits),
                          linetype='dashed', color="grey31")

    labels <- data.frame(x=max(autosomes$ID),
                         y=max(autosomes$NormHits) * 1.010,
                         label='Female')
    pp <- pp + geom_text(data=labels, aes(x=x, y=y, label=label),
                         vjust=0, hjust=1, color="grey31")

    pp <- pp + geom_hline(yintercept=0.5 * range(autosomes$NormHits),
                          linetype='dashed', color="grey")

    labels <- data.frame(x=max(autosomes$ID),
                         y=max(autosomes$NormHits) * 0.505,
                         label='Male?', color="grey")
    pp <- pp + geom_text(data=labels, aes(x=x, y=y, label=label),
                         vjust=0, hjust=1, color="grey")


    pp <- pp + geom_point(data=autosomes, aes(x=ID, y=NormHits))

    pp <- pp + ylab("Estimated ploidy")
    pp <- pp + xlab("Chromosome")
    pp <- pp + theme_bw()
    pp <- pp + theme(axis.ticks.x=element_blank(),
                     panel.border=element_blank())

    pp <- pp + scale_x_continuous(limits=range(autosomes$ID),
                                  breaks=seq(1, max(autosomes$ID) + 10, 10))


    sex <- tbl[tbl$ID == 'X' | tbl$ID == 'Y', , drop=FALSE]
    sex <- sex[order(sex$ID), , drop=FALSE]

    if (nrow(sex) > 0) {
        id_range <- range(autosomes$ID)
        step <- (id_range[2] - id_range[1]) / (nrow(sex) + 1)
        sex$x <- id_range[1] + step * 1:nrow(sex)

        sex$NormHits <- 2 * sex$RelHits / mean(autosomes$RelHits)

        pp <- pp + geom_point(data=sex, shape=sex$ID, color="red", size=5,
                          aes(x=x, y=NormHits))

        ymin <- min(0.500, 0.95 * min(autosomes$NormHits, sex$NormHits))
        ymax <- max(2.500, 1.05 * max(autosomes$NormHits, sex$NormHits))

        pp <- pp + scale_y_continuous(breaks=seq(1, 2),
                                      limits=c(ymin, ymax))
    }

    return(pp)
}


input_file <- args[1]
output_prefix <- args[2]

pdf(paste(output_prefix, ".pdf", sep=""), width=5, height=5)
plot_coverage(input_file)
dev.off()

# bitmap is preferred, since it works in a headless environment
bitmap(paste(output_prefix, ".png", sep=""), height=5, width=5, res=96, taa=4, gaa=4)
plot_coverage(input_file)
dev.off()
