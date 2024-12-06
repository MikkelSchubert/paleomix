#!/bin/env Rscript
#
# Copyright (c) 2009-2020, Broad Institute, Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name Broad Institute, Inc. nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
library(tools)

args <- commandArgs(TRUE)
if (length(args) < 2) {
    cat("Usage: tranches.r <input> <output_prefix> [TITV] [sensitivity] [suppressLegend]")
    quit(status=1)
}

tranchesFile = args[1]
outputPrefix = args[2]
targetTITV = as.numeric(args[3])
targetSensitivity = as.numeric(args[4])
suppressLegend = !is.na(args[5])

if (is.na(targetTITV)) {
    targetTITV = 2
}

# -----------------------------------------------------------------------------------------------
# Useful general routines
# -----------------------------------------------------------------------------------------------

MIN_FP_RATE = 0.001  # 1 / 1000 is min error rate

titvFPEst <- function(titvExpected, titvObserved) {
    max(min(1 - (titvObserved - 0.5)/(titvExpected - 0.5), 1), MIN_FP_RATE)
}

titvFPEstV <- function(titvExpected, titvs) {
    sapply(titvs, function(x) titvFPEst(titvExpected, x))
}

nTPFP <- function(nVariants, FDR) {
    return(list(TP = nVariants * (1 - FDR/100), FP = nVariants * (FDR/100)))
}

leftShift <- function(x, leftValue = 0) {
    r = rep(leftValue, length(x))
    for (i in 1:(length(x) - 1)) {
        # print(list(i=i))
        r[i] = x[i + 1]
    }
    r
}

# -----------------------------------------------------------------------------------------------
# Tranches plot
# -----------------------------------------------------------------------------------------------
data2 = read.table(tranchesFile, sep = ",", head = T)
# data2 = data2[order(data2$novelTiTv, decreasing = F), ]
# data2 = data2[order(data2$FDRtranche, decreasing=T),]
data2 = data2[order(data2$truthSensitivity, decreasing = T), ]
cols = c("cornflowerblue", "cornflowerblue", "darkorange", "darkorange")
density = c(20, -1, -1, 20)
outfile = paste(outputPrefix, ".pdf", sep = "")
pdf(outfile, height = 5, width = 8)
par(mfrow = 1:2)
par(mar = c(5, 5, 4, 2) + 0.1)
novelTiTv = c(data2$novelTITV, data2$novelTiTv)
alpha = 1 - titvFPEstV(targetTITV, novelTiTv)
# print(alpha)

numGood = round(alpha * data2$numNovel)

# numGood = round(data2$numNovel * (1-data2$targetTruthSensitivity/100))
numBad = data2$numNovel - numGood

numPrevGood = leftShift(numGood, 0)
numNewGood = numGood - numPrevGood
numPrevBad = leftShift(numBad, 0)
numNewBad = numBad - numPrevBad

d = matrix(c(numPrevGood, numNewGood, numNewBad, numPrevBad), 4, byrow = TRUE)
# print(d)
barplot(d/1000, horiz = TRUE, col = cols, space = 0.2, xlab = "Number of Novel Variants (1000s)",
    density = density, cex.axis = 1.25, cex.lab = 1.25)  # , xlim=c(250000,350000))
# abline(v= d[2,dim(d)[2]], lty=2) abline(v= d[1,3], lty=2)
if (!suppressLegend) {
    legend("topright", c("Cumulative TPs", "Tranch-specific TPs", "Tranch-specific FPs",
        "Cumulative FPs"), fill = cols, density = density, bg = "white", cex = 0.8)
}

mtext("Ti/Tv", 2, line = 2.25, at = length(data2$targetTruthSensitivity) * 1.2, las = 1,
    cex = 1)
mtext("truth", 2, line = 0, at = length(data2$targetTruthSensitivity) * 1.2, las = 1,
    cex = 1)
axis(2, line = -1, at = 0.7 + (0:(length(data2$targetTruthSensitivity) - 1)) * 1.2,
    tick = FALSE, labels = data2$targetTruthSensitivity, las = 1, cex.axis = 1)
axis(2, line = 1, at = 0.7 + (0:(length(data2$targetTruthSensitivity) - 1)) * 1.2,
    tick = FALSE, labels = round(novelTiTv, 3), las = 1, cex.axis = 1)

# plot sensitivity vs. specificity
sensitivity = data2$truthSensitivity
if (!is.null(sensitivity)) {
    # specificity = titvFPEstV(targetTITV, novelTiTv)
    specificity = novelTiTv
    plot(sensitivity, specificity, type = "b", col = "cornflowerblue", xlab = "Tranche truth sensitivity",
        ylab = "Specificity (Novel Ti/Tv ratio)")
    abline(h = targetTITV, lty = 2)
    abline(v = targetSensitivity, lty = 2)
    # text(max(sensitivity), targetTITV-0.05, labels='Expected novel Ti/Tv', pos=2)
}

dev.off()

if (exists("compactPDF")) {
    compactPDF(outfile)
}

res <- as.data.frame(t(d))
names(res) <- c("CumTPs", "TranchSpecificTPs", "TranchSpecificFPs", "CumFPs")
rownames(res) <- colnames(d)
res$totalNovel <- res$CumTPs + res$CumFPs
res$novelTiTv <- novelTiTv
res$tranches <- sensitivity

outfile = paste(outputPrefix, ".txt", sep = "")
write.table(res, file = outfile, col = T, row = F, qu = F, sep = "\t")
