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
require(ggplot2)
require(reshape2)


plot.admixture <- function(input_file, sample_names)
{
    samples <- read.table(sample_names, as.is=TRUE, comment.char="", header=TRUE)
    Q <- read.table(input_file)
    Q <- t(as.matrix(Q))
    colnames(Q) <- samples$Name

    # Order by name, and then move Sample to the right (Group is '-')
    Q <- Q[, order(samples$Group, samples$Name)]
    Q <- cbind(Q[, -1], Q[, 1, drop=FALSE])
    color <- samples$Color[order(samples$Group, samples$Name)]
    color <- c(color[-1], color[1])

    pp <- ggplot(melt(Q), aes(x=Var2, y=value, fill=Var1))
    pp <- pp + geom_bar(stat="identity")

    pp <- pp + xlab(NULL)
    pp <- pp + ylab("Ancestry")

    pp <- pp + theme_minimal()
    pp <- pp + theme(legend.position="None",
                     axis.text.y=element_blank(),
                     panel.grid.minor.y=element_blank(),
                     panel.grid.major.y=element_blank(),
                     axis.text.x=element_text(angle=25,
                                              hjust=1,
                                              size=12,
                                              color=color))

    return(pp)
}



args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    cat("Usage: admixture.R <input_file> <sample_names> <output_prefix>\n", file=stderr())
    quit(status=1)
}

input_file <- args[1]
sample_names <- args[2]
output_prefix <- args[3]

pdf(paste(output_prefix, ".pdf", sep=""))
plot.admixture(input_file, sample_names)
dev.off()

# bitmap is preferred, since it works in a headless environment
bitmap(paste(output_prefix, ".png", sep=""), height=6, width=6, res=96, taa=4, gaa=4)
plot.admixture(input_file, sample_names)
dev.off()
