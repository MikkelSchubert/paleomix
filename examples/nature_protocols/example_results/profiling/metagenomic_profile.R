############################################################
## Part 46:
input <- read.csv("Potato_merged.csv", header=TRUE, row.names=1, sep="\t")

dim(input)

abundances <- input[row.names(input)[grep("g__[\\w_]*$", row.names(input), perl=TRUE)],]
row.names(abundances) <- sapply(strsplit(row.names(abundances), "g__"), `[`, 2)
write.csv(abundances, "Relative_abundances_genus.csv")


############################################################
## Part 47:

# Load required packages
library(permute)
library(vegan)
library(MASS)
# Calculate the number of genera identified for each profile
taxon_table <- specnumber(t(abundances))
# Export a summary table
write.table(taxon_table, "Taxon_count_genus.txt", sep="\t", col.names=FALSE)
# Calculate the Shannon diversity index for each profile at the genus level
diversity_table <- diversity(t(abundances), index="shannon")
# Export a summary table
write.table(diversity_table, "Diversity_genus.txt", sep="\t", col.names=FALSE)
# Calculate Bray-Curtis distances among profiles at the genus level
distances <- vegdist(t(abundances), method="bray")
# Export a summary table
write.matrix(as.matrix(distances), "Distances_genus.txt", sep="\t")


############################################################
## Part 48:

# Load required packages
library(gplots)
# Save heatmap as .pdf file
pdf("Heatmap_genus.pdf")
# Draw the heatmap (see below)
heatmap.2(as.matrix(abundances), trace="none", col=colorpanel(50, low="gray91", mid="darkblue", high="red"), dendrogram="both", ColSideColors=rep(c("red","green","blue"), times=c(5,1,3)), margins=c(12, 13))
# Close the .pdf file
dev.off()


############################################################
## Part 49:

# Perform the Principal Coordinate Analysis
library(ape)
pcoa <- pcoa(distances)
# Save the PCA plot as a pdf file.
pdf("PCOA_genus.pdf")


############################################################
## Part 50:

# Plot the first two dimensions
plot(pcoa$vectors[,1], pcoa$vectors[,2], pch=16, cex=1.5, cex.axis=0.9, font.lab=2, font.axis=2, xlab="Dimension 1", ylab="Dimension 2", col=rep(c("red","green","blue"), times=c(5, 1, 3)))
# Add profile name labels, colored using the same color scheme as above.
text(x=pcoa$vectors[,1], y=pcoa$vectors[,2], labels=row.names(pcoa$vectors), font=2, cex=0.8, pos=3, col=rep(c("red","green","blue"), times=c(5, 1, 3)))
# Add a legend showing the correspondence between profiles and samples.
legend('topright',  legend=c("M-0182896", "Pi1845A", "Pi1889"), pch=16, col=c('red', 'green', 'blue'), bty='n', cex=.75)
dev.off()

# Perform Principal Component Analysis
pca <- prcomp(t(abundances), scale.=T)
# Calculate the percentage of the variance accounted by principal components 1 and 2
PC1 <- round (100 * (summary(pca)$importance[2,1]), 1)
PC2 <- round (100 * (summary(pca)$importance[2,2]), 1)
# Plot PCA scores for principal component 1 and 2
pdf("PCA_genus.pdf")

# Plot the first two principal components
plot(pca$x[,1], pca$x[,2], pch=16, cex=1.5, cex.axis=0.9, font.lab=2, font.axis=2, xlab=paste("Principal component 1 - ", PC1, "% variance", sep=""), ylab=paste("Principal component 2 - ", PC2, "% variance", sep=""), col=rep(c("red","green","blue"), times=c(5, 1, 3)))
# Add profile name labels, colored using the same color scheme as above.
text(pca$x[,1], pca$x[,2], colnames(abundances), font=2, cex=0.8, pos=3, col=rep(c("red", "green", "blue"), times=c(5, 1, 3)))
# Add a legend showing the correspondence between profiles and samples.
legend('topright',  legend=c("M-0182896", "Pi1845A", "Pi1889"), pch=16, col=c('red', 'green', 'blue'), bty='n', cex=.75)

# Plot PCA loadings and their labels
vectors.x <- (pca$rotation[,1]) * 8
vectors.y <- (pca$rotation[,2]) * 8
points(cbind(vectors.x, vectors.y), col="grey50", type="n")
text(cbind(vectors.x, vectors.y), rownames(cbind(vectors.x, vectors.y)), cex=0.6, font=2, pos=3, col="grey50")
for (v in 1:length(vectors.x)) {
    segments(0,0, vectors.x[v], vectors.y[v], col="grey50", lty=3, lwd=2.5)
}
dev.off()


############################################################
## Part 51:

library(pvclust)
clustering <- pvclust(abundances, method.dist="manhattan", method.hclust="average", nboot=1000)
pdf("Clustering_genus.pdf")
plot(clustering)
dev.off()

clustering


############################################################
## Part 53:

# Import and format class-level abundance data
input2 <- read.csv("Potato_merged_profiles.csv", header=TRUE, row.names=1, sep="\t")
abundances2 <- input2[row.names(input2)[grep("c__[\\w_]*$", row.names(input2), perl=TRUE)],]
row.names(abundances2) <- sapply(strsplit(row.names(abundances2), "c__"), `[`, 2)
data_table <- data.frame(samples=rep(colnames(abundances2), each=dim(abundances2)[1]), taxa=rep(row.names(abundances2), dim(abundances2)[2]), datavect=unlist(abundances2))

library(ggplot2)
library(grid)
# Draw the stacked bar plot
ggplot(data_table, aes(x=samples)) + geom_bar(aes(weight=datavect, fill=taxa), position='fill') + scale_y_continuous("", breaks=NULL) + scale_fill_manual(values=rainbow(dim(abundances2)[1])) + theme(axis.text.x=element_text(angle=90, hjust = 0, color="black"), legend.text=element_text(size = 8)) + theme(legend.key.size=unit(0.6, "lines"))
ggsave("Barplot_class.pdf")
