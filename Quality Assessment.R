##################################################
##################################################

# Have just the Cel files in one folder on your desktop

cel.folder.name <- "Microarray_CEL files"
# for example: cel.folder.name <- "Cel Files"
# must type the name exactly as it appears (capitalization, spaces, etc)

##################################################
##################################################

library(oligo)

dir.create(paste("~/Desktop", "PUFA Microarray Analysis", sep="/"))
main.directory <- paste("~/Desktop", "PUFA Microarray Analysis", sep="/")
dir.create(paste(main.directory, "All 32 Arrays", sep="/"))
subdir.all <- paste(main.directory, "All 32 Arrays/", sep="/")

setwd(paste("~/Desktop", cel.folder.name, sep="/"))

array.list <- list.files()
raw <- read.celfiles(array.list)


##### Name Files #####

name.frame <- data.frame(LC=paste("LC", 1:8, sep=""), OC=paste("OC", 1:8, sep=""), ALA=paste("ALA", 1:8, sep=""), LA=paste("LA", 1:8, sep=""))

t.name.frame <- t(name.frame)

condition.names <- as.vector(c(t.name.frame[,1], t.name.frame[,2], t.name.frame[,3], t.name.frame[,4], t.name.frame[,5], t.name.frame[,6], t.name.frame[,7], t.name.frame[,8]))

sampleNames(raw) <- c(condition.names)
# applies condition names and replicate numbers to each array


##### Quality Assessment #####

plot.colors <- c(rep(c("red", "green", "cyan", "purple"), 8))
group.colors <- c(rep("red", 8), rep("green", 8), rep("cyan", 8), rep("purple", 8))
range.colors <- rainbow(32)

group <- c(seq(from=1, to=32, by=4), seq(from=2, to=32, by=4), seq(from=3, to=32, by=4), seq(from=4, to=32, by=4))
# subsetting vector to group treatment groups together

dir.create(paste(subdir.all, "QA Images of Raw Data", sep="/"))
subdir.allQA <- paste(subdir.all, "QA Images of Raw Data/", sep="/")


### Boxplots ###

pdf(file=paste(subdir.allQA, "Raw Boxplot.pdf", sep=""), width=20, height=7)
boxplot(raw, range=1.5, col=plot.colors, xlab="Array", ylab="Log Probe Intensity", main="Raw Log Probe Intensity")
dev.off()

pdf(file=paste(subdir.allQA, "Raw Boxplot Grouped.pdf", sep=""), width=20, height=7)
boxplot(raw[,group], range=1.5, col=group.colors, xlab="Array", ylab="Log Probe Intensity", main="Raw Log Probe Intensity")
dev.off()


### Density Plot ###

pdf(file=paste(subdir.allQA, "Raw Density Estimation Plot.pdf", sep=""), width=10, height=7)
hist(raw, col=range.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=range.colors, lty=1)
dev.off()

pdf(file=paste(subdir.allQA, "Raw Density Estimation Plot Grouped.pdf", sep=""), width=10, height=7)
hist(raw[,group], col=group.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names[group]), col=group.colors, lty=1)
dev.off()


### Principal Component Analysis Plot ###

pca.legend.colors <- c("red", "green", "cyan", "purple")
pca.conditions <- c("LC", "OC", "ALA", "LA")
pca.numbers <- c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 4), rep(7, 4), rep(8, 4))

raw.expression.matrix <- exprs(raw)
transposed.raw.expression.matrix <- t(raw.expression.matrix)

pca.values.raw <- prcomp(transposed.raw.expression.matrix)

normal<-par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=F)

pdf(file=paste(subdir.allQA, "Raw PCA Plot.pdf", sep=""))
par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
# adds space on the side of the graph for the legend

plot(pca.values.raw$x, col=plot.colors, pch=20, main="Raw PCA Plot")
text(pca.values.raw$x, pos=3, offset=0.2, labels=pca.numbers, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
dev.off()

par(normal)
# returns the graph parameters to normal

pca.summary.raw <- summary(pca.values.raw)
proportion.variance.raw <- pca.summary.raw$importance[2:3,1:5]
# gets the values corresponding to the proportion of variance and cumulative variance for the first five principal components

pdf(file=paste(subdir.allQA, "Raw Proportion of Variance PCA.pdf", sep=""), width=7, height=7)
barplot(proportion.variance.raw, beside=T, col=c("black","gray"), main="Raw Proportion of Variance of Principal Components", xlab="Principal Components", ylab="Proportion")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
dev.off()


### Hierarchical Clustering Dendogram ###

transposed <- t(raw.expression.matrix)
distance <- dist(transposed)
sample.clusters <- hclust(distance)

pdf(file=paste(subdir.allQA, "Raw Hierarchical Cluster Dendogram.pdf", sep=""))
plot(sample.clusters, main="Raw Hierarchical Cluster Dendogram", xlab="Samples", sub="")
dev.off()


### Image Plots ###

dir.create(paste(subdir.allQA, "Image Plots", sep="/"))
subdir.allQA.im <- paste(subdir.allQA, "Image Plots/", sep="/")

for (i in 1:32){
	png(file=paste(subdir.allQA.im, paste(condition.names[i], "Image Plot.png", sep=" "), sep=""))
	image(raw[,i], xlab="Probe Intensity")
	dev.off()
}


### Log-Transformed Image Plots ###

dir.create(paste(subdir.allQA, "Log Image Plots", sep="/"))
subdir.allQA.ltim <- paste(subdir.allQA, "Log Image Plots/", sep="/")

for (i in 1:32){
	png(file=paste(subdir.allQA.ltim, paste(condition.names[i], "Log Image Plot.png", sep=" "), sep=""))
	image(raw[,i], xlab="Log Probe Intensity", transfo=function(x)x)
	dev.off()
}


### Probe Level Model Fit Quality Assessment ###

raw.plm <- fitProbeLevelModel(raw)

dir.create(paste(subdir.allQA, "PLM Plots", sep="/"))
subdir.allQA.plm <- paste(subdir.allQA, "PLM Plots/", sep="/")


# NUSE Plots #

pdf(file=paste(subdir.allQA.plm, "NUSE plot.pdf", sep=""), width=20, height=7)
NUSE(raw.plm, xlab="Array", main="Normalized Unscaled Standard Errors", col=plot.colors)
dev.off()


# RLE Plots #

pdf(file=paste(subdir.allQA.plm, "RLE plot.pdf", sep=""), width=20, height=7)
RLE(raw.plm, xlab="Array", main="Relative Log Expression", col=plot.colors)
dev.off()


# Weights Image Plots #

dir.create(paste(subdir.allQA.plm, "Weights Image Plots", sep="/"))
subdir.allQA.plmW <- paste(subdir.allQA.plm, "Weights Image Plots/", sep="/")

for (i in 1:32){
	png(file=paste(subdir.allQA.plmW, paste(condition.names[i], "PLM Weights Image.png", sep=" "), sep=""))
	image(raw.plm, type="weights", which=i, xlab="Weights PLM Plot")
	dev.off()
}


# Residuals Image Plots #

dir.create(paste(subdir.allQA.plm, "Residuals Image Plots", sep="/"))
subdir.allQA.plmR <- paste(subdir.allQA.plm, "Residuals Image Plots/", sep="/")

for (i in 1:32){
	png(file=paste(subdir.allQA.plmR, paste(condition.names[i], "PLM Residuals Image.png", sep=" "), sep=""))
	image(raw.plm, type="residuals", which=i, xlab="Residuals PLM Plot")
	dev.off()
}


# Signs of the Residuals Image Plots #

dir.create(paste(subdir.allQA.plm, "Signs of Residuals Image Plots", sep="/"))
subdir.allQA.plmS <- paste(subdir.allQA.plm, "Signs of Residuals Image Plots/", sep="/")

for (i in 1:32){
	png(file=paste(subdir.allQA.plmS, paste(condition.names[i], "PLM Signs Image.png", sep=" "), sep=""))
	image(raw.plm, type="sign.residuals", which=i, xlab="Signs of the Residuals PLM Plot")
	dev.off()
}


##### Preprossess the Raw Data #####

normalized <- rma(raw)

dir.create(paste(subdir.all, "Preprocessed Data", sep="/"))
subdir.all.preproc <- paste(subdir.all, "Preprocessed Data/", sep="/")
dir.create(paste(subdir.all.preproc, "Images", sep="/"))
subdir.all.preproc.im <- paste(subdir.all.preproc, "Images/", sep="/")


### Boxplot ###

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Boxplot.pdf", sep=""), width=20, height=7)
boxplot(normalized, range=1.5, col=plot.colors, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity")
dev.off()

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Boxplot Grouped.pdf", sep=""), width=20, height=7)
boxplot(normalized[,group], range=1.5, col=group.colors, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity")
dev.off()


### Density Plot ###

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Density Estimation Plot.pdf", sep=""), width=10, height=7)
hist(normalized, col=range.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Preprocessed Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=range.colors, lty=1)
dev.off()

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Density Estimation Plot Grouped.pdf", sep=""), width=10, height=7)
hist(normalized[,group], col=group.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Preprocessed Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names[group]), col=group.colors, lty=1)
dev.off()


### Principal Component Analysis ###

preprocessed.expression.matrix <- exprs(normalized)
transposed.preprocessed.expression.matrix <- t(preprocessed.expression.matrix)

pca.values.preprocessed <- prcomp(transposed.preprocessed.expression.matrix)

pdf(file=paste(subdir.all.preproc.im, "Preprocessed PCA Plot.pdf", sep=""))
par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.preprocessed$x, col=plot.colors, pch=20, main="Preprocessed PCA Plot")
text(pca.values.preprocessed$x, pos=3, offset=0.2, labels=pca.numbers, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
par(normal)
dev.off()

pca.summary.preprocessed <- summary(pca.values.preprocessed)
proportion.variance.preprocessed <- pca.summary.preprocessed$importance[2:3,1:5]

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Proportion of Variance PCA.pdf", sep=""))
barplot(proportion.variance.preprocessed, beside=T, col=c("black","gray"), main="Preprocessed Proportion of Variance of Principal Components", xlab="Principal Components", ylab="Proportion")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
dev.off()


### Hierarchical Clustering Dendogram ###

preprocessed.transposed <- t(preprocessed.expression.matrix)
preprocessed.distance <- dist(preprocessed.transposed)
preprocessed.sample.clusters <- hclust(preprocessed.distance)

pdf(file=paste(subdir.all.preproc.im, "Preprocessed Hierarchical Cluster Dendogram.pdf", sep=""))
plot(preprocessed.sample.clusters, main="Preprocessed Hierarchical Cluster Dendogram", xlab="Samples")
dev.off()
