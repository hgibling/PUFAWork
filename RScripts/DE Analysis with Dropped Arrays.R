##################################################
##################################################

# If after performing Quality Assessment it is decided some arrays should be dropped: 

# Have the annotation file (Annotations for Rat Gene 21st.csv) on your desktop

arrays.to.be.dropped <- c(18)
# for example: arrays.to.be.dropped<-c(5,6)

# Input remaining arrays for each condition and NA for dropped arrays

LC.remaining <- c(1:8)
OC.remaining <- c(1:4, 6:8)
ALA.remaining <- c(1:8)
LA.remaining <- c(1:8)
# for example: LC.remaining <- c(1, NA 3:8)
# if array 5 has been dropped, the remaining LC samples are LC1, LC3 through LC8

##################################################
##################################################


dir.create(paste(main.directory, "Data With Dropped Arrays (OC5)", sep="/"))
subdir.drop <- paste(main.directory, "Data With Dropped Arrays (OC5)/", sep="/")

raw.dropped <- raw[,-c(arrays.to.be.dropped)]

condition.names.dropped <- condition.names[-arrays.to.be.dropped]

sampleNames(raw.dropped) <- c(condition.names.dropped)
# applies condition names and replicate numbers to each array


##### Quality Assessment with Dropped Arrays #####

plot.colors.dropped <- plot.colors[-arrays.to.be.dropped]
group.colors.dropped <- group.colors[-arrays.to.be.dropped]
range.colors.dropped <- rainbow(as.integer(dim(raw.dropped)[2]))

drop.position <- which(group %in% arrays.to.be.dropped)
group.dropped.fix <- group[-drop.position]
group.dropped <- group.dropped.fix

for (i in 1:length(group.dropped.fix)){
	if (group.dropped.fix[i] > arrays.to.be.dropped[1]){
		group.dropped[i] <- (group.dropped.fix[i] - 1)
	} else {
		group.dropped[i] <- group.dropped.fix[i]
	}
}
# only works when one array has been dropped... need to fix for multiple dropped arrays

dir.create(paste(subdir.drop, "QA Images of Raw Data", sep="/"))
subdir.dropQA <- paste(subdir.drop, "QA Images of Raw Data/", sep="/")


### Boxplots ###

pdf(file=paste(subdir.dropQA, "Raw Boxplot with Dropped Arrays.pdf", sep=""), width=20, height=7)
boxplot(raw.dropped, range=1.5, col=plot.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main="Raw Log Probe Intensity with Dropped Arrays")
dev.off()

pdf(file=paste(subdir.dropQA, "Raw Boxplot with Dropped Arrays Grouped.pdf", sep=""), width=20, height=7)
boxplot(raw.dropped[,group.dropped], range=1.5, col=group.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main="Raw Log Probe Intensity with Dropped Arrays")
dev.off()


### Density Plot ###

hist(raw.dropped, col=range.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation with Dropped Arrays")
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped), col=plot.colors.dropped, lty=1)
quartz.save(paste(subdir.dropQA, "Raw Density Estimation Plot with Dropped Arrays.pdf", sep=""), type="pdf", width=10, height=7)

pdf(file=paste(subdir.dropQA, "Raw Density Estimation Plot with Dropped Arrays.pdf", sep=""), width=10, height=7)
hist(raw.dropped, col=range.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation with Dropped Arrays")
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped), col=range.colors, lty=1)
dev.off()

pdf(file=paste(subdir.dropQA, "Raw Density Estimation Plot  with Dropped Arrays Grouped.pdf", sep=""), width=10, height=7)
hist(raw.dropped[,group.dropped], col=group.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation with Dropped Arrays")
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped[group.dropped]), col=group.colors.dropped, lty=1)
dev.off()



### Principal Component Analysis Plot ###

pca.numbers.dropped <- pca.numbers[-arrays.to.be.dropped]

raw.expression.matrix.dropped <- exprs(raw.dropped)
transposed.raw.expression.matrix.dropped <- t(raw.expression.matrix.dropped)

pca.values.raw.dropped <- prcomp(transposed.raw.expression.matrix.dropped)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.raw.dropped$x, col=pca.colors.dropped, pch=20, main="Raw PCA Plot with Dropped Arrays")
text(pca.values.raw.dropped$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.dropQA, "Raw PCA Plot with Dropped Arrays.pdf", sep=""), type="pdf")
par(normal)

pca.summary.raw.dropped <- summary(pca.values.raw.dropped)
proportion.variance.raw.dropped <- pca.summary.raw.dropped$importance[2:3,1:5]

barplot(proportion.variance.raw.dropped, beside=T, col=c("black","gray"), main="Raw Proportion of Variance of Principal Components \n with Dropped Arrays", xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.dropQA, "Raw Proportion of Variance PCA with Dropped Arrays.pdf", sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

transposed.dropped <- t(raw.expression.matrix.dropped)
distance.dropped <- dist(transposed.dropped)
sample.clusters.dropped <- hclust(distance.dropped)

plot(sample.clusters.dropped, main="Raw Hierarchical Cluster Dendogram \n with Dropped Arrays", xlab="Samples", sub="")
quartz.save(paste(subdir.dropQA, "Raw Hierarchical Cluster Dendogram with Dropped Arrays.pdf", sep=""), type="pdf")


### PLM Fit ###

raw.dropped.plm <- fitProbeLevelModel(raw.dropped)

dir.create(paste(subdir.dropQA, "PLM Plots", sep="/"))
subdir.dropQA.plm <- paste(subdir.dropQA, "PLM Plots/", sep="/")


# NUSE Plots #

NUSE(raw.dropped.plm, xlab="Array", main="Normalized Unscaled Standard Errors with Dropped Arrays")
quartz.save(paste(subdir.dropQA.plm, "NUSE plot with Dropped Arrays.pdf", sep=""), type="pdf", width=15, height=7)


# RLE Plots #

RLE(raw.dropped.plm, xlab="Array", main="Relative Log Expression with Dropped Arrays")
quartz.save(paste(subdir.dropQA.plm, "RLE plot with Dropped Arrays.pdf", sep=""), type="pdf", width=15, height=7)


##### Preprocess with Dropped Arrays #####

normalized.dropped <- rma(raw.dropped)

dir.create(paste(subdir.drop, "Preprocessed Data", sep="/"))
subdir.drop.preproc <- paste(subdir.drop, "Preprocessed Data/", sep="/")
dir.create(paste(subdir.drop.preproc, "Images", sep="/"))
subdir.drop.preproc.im <- paste(subdir.drop.preproc, "Images/", sep="/")


### Boxplot ###

boxplot(normalized.dropped, range=1.5, col=range.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity with Dropped Arrays")
quartz.save(paste(subdir.drop.preproc.im, "Preprossed Boxplot with Dropped Arrays.pdf", sep=""), type="pdf", width=15, height=7)


### Density Plot ###

hist(normalized.dropped, col=range.colors.dropped, lty=1, xlab="Log Intensity", ylab="Density", main="Preprocessed Density Estimation with Dropped Arrays")
legend("topright", inset=0.01, cex=0.75, c(condition.names.dropped), col=plot.colors.dropped, lty=1)
quartz.save(paste(subdir.drop.preproc.im, "Preprocessed Density Estimation Plot with Dropped Arrays.pdf", sep=""), type="pdf", width=10, height=7)


### Principal Component Analysis ###

preprocessed.expression.matrix.dropped <- exprs(normalized.dropped)
transposed.preprocessed.expression.matrix.dropped <- t(preprocessed.expression.matrix.dropped)

pca.values.preprocessed.dropped <- prcomp(transposed.preprocessed.expression.matrix.dropped)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.values.preprocessed.dropped$x, col=pca.colors.dropped, pch=20, main="Preprocessed PCA Plot with Dropped Arrays")
text(pca.values.preprocessed.dropped$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subdir.drop.preproc.im, "Preprocessed PCA Plot with Dropped Arrays.pdf", sep=""), type="pdf")
par(normal)

pca.summary.preprocessed.dropped <- summary(pca.values.preprocessed.dropped)
proportion.variance.preprocessed.dropped <- pca.summary.preprocessed.dropped$importance[2:3,1:5]

barplot(proportion.variance.preprocessed.dropped, beside=T, col=c("black","gray"), main="Preprocessed Proportion of Variance of Principal Components \n with Dropped Arrays", xlab="Principal Components", ylab="Percentage")
legend("topleft", inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subdir.drop.preproc.im, "Preprocessed Proportion of Variance PCA with Dropped Arrays.pdf", sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

preprocessed.transposed.dropped <- t(preprocessed.expression.matrix.dropped)
preprocessed.distance.dropped <- dist(preprocessed.transposed.dropped)
preprocessed.sample.clusters.dropped <- hclust(preprocessed.distance.dropped)

plot(preprocessed.sample.clusters.dropped, main="Preprocessed Hierarchical Cluster Dendogram \n with Dropped Arrays", xlab="Samples", sub="")
quartz.save(paste(subdir.drop.preproc.im, "Preprocessed Hierarchical Cluster Dendogram with Dropped Arrays.pdf", sep=""), type="pdf")


##### Filter Genes to Remove Those with No Annotation #####

dir.create(paste(subdir.drop.preproc, "Filtered Genes Data", sep="/"))
subdir.drop.preproc.filt <- paste(subdir.drop.preproc, "Filtered Genes Data/", sep="/")

annotation.file <- read.csv("~/Desktop/Annotations for Rat Gene 21st.csv", header=TRUE)
annotation.matrix <- as.matrix(annotation.file)
annotations <- annotation.matrix[order(annotation.file[,1], annotation.file[,2]),]
# creates a matrix where the genes are ordered by their transcript cluster

annotated.values <- cbind(annotations, preprocessed.expression.matrix.dropped)
# combines gene annotations with expression values, which are already ordered according to transcript cluster
annotated.gene.values.all <- annotated.values[,-1]
# removes transcript cluster column
annotated.gene.values.no.NA <- na.omit(annotated.gene.values.all)
# removes all genes with no annotation

rownames(annotated.gene.values.no.NA) <- annotated.gene.values.no.NA[,1]
annotated.gene.values <- annotated.gene.values.no.NA[,-1]
# removes column of gene IDs so all columns are expression values


##### Average Expression Values for Duplicated Genes #####

genes <- as.numeric(rownames(annotated.gene.values))
# warning is ok

duplicates <- data.frame(Gene=genes, annotated.gene.values)
# warning is ok

remove.duplicates <- aggregate(duplicates, by=list(duplicates[,1]), FUN=mean)
# for genes that are duplicated, for each column, the values are averaged

rownames(remove.duplicates)<-remove.duplicates[,2]

no.duplicates <- remove.duplicates[,-c(1:2)]


### Create Reference List of All Annotated Genes for FunNet ##

reference.list <- as.numeric(rownames(no.duplicates))

write.table(reference.list, paste(subdir.drop.preproc.filt, "Reference Gene List.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Images for Filtered Genes Data #####

dir.create(paste(subdir.drop.preproc.filt, "Images", sep="/"))
subsubdir.images <- paste(subdir.drop.preproc.filt, "Images/", sep="/")


### Boxplot ###

boxplot(no.duplicates, range=1.5, col=plot.colors.dropped, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity with Annotated Genes")
quartz.save(paste(subsubdir.images, "Preprossed Boxplot with Annotated Genes.pdf", sep=""), type="pdf", width=15, height=7)


### PCA Plot of Filtered Genes ###

transposed.annotated.gene.values <- t(no.duplicates)
pca.annotated.genes <- prcomp(transposed.annotated.gene.values)

par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)
plot(pca.annotated.genes$x, col=pca.colors.dropped, pch=20, main="Preprocessed PCA Plot with Annotated Genes")
text(pca.annotated.genes$x, pos=3, offset=0.2, labels=pca.numbers.dropped, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
quartz.save(paste(subsubdir.images, "Preprocessed PCA Plot with Annotated Genes", sep=""), type="pdf")
par(normal)

pca.summary <- summary(pca.annotated.genes)
proportion.variance <- pca.summary$importance[2:3,1:5]

barplot(proportion.variance,beside=T, col=c("black","gray"), main="Preprocessed Proportion of Variance of Principal Components \n with Annotated Genes", xlab="Principal Components", ylab="Percentage")
legend("topleft",inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
quartz.save(paste(subsubdir.images, "Preprocessed Proportion of Variance PCA with Annotated Genes.pdf", sep=""), type="pdf", width=7, height=7)


### Hierarchical Clustering Dendogram ###

annotated.transposed <- t(no.duplicates)
annotated.distance <- dist(annotated.transposed)
annotated.sample.clusters <- hclust(annotated.distance)

plot(annotated.sample.clusters, main="Preprocessed Hierarchical Cluster Dendogram \n with Annotated Genes", xlab="Samples", sub="")
quartz.save(paste(subsubdir.images, "Preprocessed Hierarchical Cluster Dendogram with Annotated Genes.pdf", sep=""), type="pdf")


##### Rearrange Samples so Replicates are Together #####

LC.positions <- seq(1, 32, 4)
OC.positions <- LC.positions+1
ALA.positions <- LC.positions+2
LA.positions <- LC.positions+3

for (i in 1:length(LC.positions)){
	if (LC.positions[i]==arrays.to.be.dropped){
		LC.positions[i] <- 0
	}
	if (LC.positions[i]>arrays.to.be.dropped){
		LC.positions[i] <- LC.positions[i]-1
	}
}

for (i in 1:length(OC.positions)){
	if (OC.positions[i]==arrays.to.be.dropped){
		OC.positions[i] <- 0
	}
	if (OC.positions[i]>arrays.to.be.dropped){
		OC.positions[i] <- OC.positions[i]-1
	}
}

for (i in 1:length(ALA.positions)){
	if (ALA.positions[i]==arrays.to.be.dropped){
		ALA.positions[i] <- 0
	}
	if (ALA.positions[i]>arrays.to.be.dropped){
		ALA.positions[i] <- ALA.positions[i]-1
	}
}

for (i in 1:length(LA.positions)){
	if (LA.positions[i]==arrays.to.be.dropped){
		LA.positions[i] <- 0
	}
	if (LA.positions[i]>arrays.to.be.dropped){
		LA.positions[i] <- LA.positions[i]-1
	}
}

LC.arrays <- no.duplicates[LC.positions]
OC.arrays <- no.duplicates[OC.positions]
ALA.arrays <- no.duplicates[ALA.positions]
LA.arrays <- no.duplicates[LA.positions]

no.duplicates <- cbind(LC.arrays, OC.arrays, ALA.arrays, LA.arrays)


##### Use ANOVA to Determine Differentially Expressed Genes #####

# code modified from Pavlidis, P. 2003 Methods 31(4):282-289

dir.create(paste(subdir.drop.preproc.filt, "Differentially Expressed Gene Lists", sep="/"))
subsubdir.DE <- paste(subdir.drop.preproc.filt, "Differentially Expressed Gene Lists/", sep="/")

remaining.conditions <- c(rep("LC", ncol(LC.arrays)), rep("OC", length(OC.arrays)), rep("ALA", length(ALA.arrays)), rep("LA", length(LA.arrays)))

conditions <- factor(remaining.conditions)
# allows you to perform one-way ANOVA with unequal sample sizes

anova.function <- function(x){
	data <- data.frame(conditions,x)
	anova(aov(x~conditions, data))
}

anova.results <- apply(no.duplicates, 1, anova.function)
# applies ANOVA function to each row (gene) of the filtered genes

pvalue.function <- function(x){
	x["Pr(>F)"][1,]
}

pvalues <- data.frame(lapply(anova.results, pvalue.function))

pvalues.table <- t(pvalues)
colnames(pvalues.table) <- "Condition P Value"


### False-Discovery Rate (FDR) Correction ###

adjusted.pvalues <- as.data.frame(p.adjust(pvalues.table[,1], method="fdr"))
colnames(adjusted.pvalues) <- "Adjusted P Value"

p.05 <- which(adjusted.pvalues[,1]<0.05)

diff.expressed.genes <- no.duplicates[p.05,]
# gets a list of genes that have a p-value of less than 0.05 -- differentially expressed

ANOVA.geneIDs <- rownames(diff.expressed.genes)

write.table(ANOVA.geneIDs, paste(subsubdir.DE, "ANOVA DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Post-Hoc Pairwise Comparisons with FDR #####

dir.create(paste(subsubdir.DE, "Pairwise DE Genes", sep="/"))
subsubdir.DE.pair <- paste(subsubdir.DE, "Pairwise DE Genes/", sep="/")

pairwise.comparisons <- c("OCvLC", "ALAvLC", "ALAvOC", "LAvLC", "LAvOC", "ALAvLA")

diff.expressed.matrix <- as.matrix(diff.expressed.genes)

ttest.results <- matrix(nrow=nrow(diff.expressed.matrix), ncol=6)
colnames(ttest.results) <- pairwise.comparisons

for (i in 1:nrow(diff.expressed.matrix)){
	p <- pairwise.t.test(diff.expressed.matrix[i,], conditions, p.adj="fdr")
	pv <- p$p.value
	ttest.results[i,1] <- pv[3,3]
	ttest.results[i,2] <- pv[2,1]
	ttest.results[i,3] <- pv[3,1]
	ttest.results[i,4] <- pv[2,2]
	ttest.results[i,5] <- pv[3,2]
	ttest.results[i,6] <- pv[1,1]
}
# Different order from regular DE Analysis script because when you manually make a vector a factor, the levels are listed alphabetically (ALA, LA, LC, OC instead of LC, OC, ALA, LA)

rownames(ttest.results) <- rownames(diff.expressed.genes)

pairwise.diff.genes <- matrix(nrow=nrow(diff.expressed.genes), ncol=6)

for (i in 1:nrow(ttest.results)){
	if (ttest.results[i,1]<=0.05){
		pairwise.diff.genes[i,1] <- (i)
	}
	if (ttest.results[i,2]<=0.05){
		pairwise.diff.genes[i,2] <- (i)
	}
	if (ttest.results[i,3]<=0.05){
		pairwise.diff.genes[i,3] <- (i)
	}
	if (ttest.results[i,4]<=0.05){
		pairwise.diff.genes[i,4] <- (i)
	}
	if (ttest.results[i,5]<=0.05){
		pairwise.diff.genes[i,5] <- (i)
	}
	if (ttest.results[i,6]<=0.05){
		pairwise.diff.genes[i,6] <- (i)
	}
}

OCvLC.info <- na.omit(pairwise.diff.genes[,1])
ALAvLC.info <- na.omit(pairwise.diff.genes[,2])
ALAvOC.info <- na.omit(pairwise.diff.genes[,3])
LAvLC.info <- na.omit(pairwise.diff.genes[,4])
LAvOC.info <- na.omit(pairwise.diff.genes[,5])
ALAvLA.info <- na.omit(pairwise.diff.genes[,6])
# gets the row numbers for the differentially expressed genes for each pairwise comparison

OCvLC.diff.values <- as.matrix(ttest.results[,1][OCvLC.info])
ALAvLC.diff.values <- as.matrix(ttest.results[,2][ALAvLC.info])
ALAvOC.diff.values <- as.matrix(ttest.results[,3][ALAvOC.info])
LAvLC.diff.values <- as.matrix(ttest.results[,4][LAvLC.info])
LAvOC.diff.values <- as.matrix(ttest.results[,5][LAvOC.info])
ALAvLA.diff.values <- as.matrix(ttest.results[,6][ALAvLA.info])
# gets the expression values for the differentially expressed genes for each pairwise comparison

OCvLC.genes <- as.numeric(rownames(OCvLC.diff.values))
ALAvLC.genes <- as.numeric(rownames(ALAvLC.diff.values))
ALAvOC.genes <- as.numeric(rownames(ALAvOC.diff.values))
LAvLC.genes <- as.numeric(rownames(LAvLC.diff.values))
LAvOC.genes <- as.numeric(rownames(LAvOC.diff.values))
ALAvLA.genes <- as.numeric(rownames(ALAvLA.diff.values))
# gets a list of just the differentially expressed gene IDs for each pairwise comparison


### Save Gene Lists to Files ###

dir.create(paste(subsubdir.DE.pair, "All DE Genes", sep="/"))
subsubdir.DE.pairA<-paste(subsubdir.DE.pair,"All DE Genes/",sep="/")

write.table(OCvLC.genes, paste(subsubdir.DE.pairA, "OCvLC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLC.genes, paste(subsubdir.DE.pairA, "ALAvLC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvOC.genes,paste(subsubdir.DE.pairA, "ALAvOC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvLC.genes,paste(subsubdir.DE.pairA, "LAvLC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvOC.genes,paste(subsubdir.DE.pairA, "LAvOC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLA.genes, paste(subsubdir.DE.pairA, "ALAvLA DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Determine Up- and Down-Regulated Genes #####

average.condition.names <- c(rep("LC", ncol(LC.arrays)), rep("OC", length(OC.arrays)), rep("ALA", length(ALA.arrays)), rep("LA", length(LA.arrays)))

condition.value.columns <- diff.expressed.matrix

colnames(condition.value.columns) <- average.condition.names
# replaces the column names representing each individual array number with general names representing the conditions, so that the expression values for each condition can be averaged (next steps)

transposed.condition.values <- t(condition.value.columns)

group.condition.values <- data.frame(Condition= average.condition.names,  transposed.condition.values)
# warning is ok

mean.condition.values <- aggregate(group.condition.values, by=list(group.condition.values[,1]), FUN=mean)
# warning is ok
# averages out the values for each condition for each gene

combined.condition.values <- t(mean.condition.values)
cvalues <- data.frame(combined.condition.values[-c(1,2),])

LC.values <- as.numeric(levels(cvalues$X1))[cvalues$X1]
OC.values <- as.numeric(levels(cvalues$X2))[cvalues$X2]
ALA.values <- as.numeric(levels(cvalues$X3))[cvalues$X3]
LA.values <- as.numeric(levels(cvalues$X4))[cvalues$X4]
# the aggragate function results in the columns being identified as factors instead of numbers; this changes each number from a level of a factor to a numeric value

average.condition.values <- data.frame(LC=LC.values, OC=OC.values, ALA=ALA.values, LA=LA.values)

gene.names <- rownames(diff.expressed.genes)

rownames(average.condition.values) <- gene.names

pairwise.differences <- data.frame(Gene=as.numeric(gene.names), OCvLC=(average.condition.values[,2]-average.condition.values[,1]), ALAvLC=(average.condition.values[,3]-average.condition.values[,1]), ALAvOC=(average.condition.values[,3]-average.condition.values[,2]), LAvLC=(average.condition.values[,4]-average.condition.values[,1]), LAvOC=(average.condition.values[,4]-average.condition.values[,2]), ALAvLA=(average.condition.values[,3]-average.condition.values[,4]))
# subtracts the expression values of the second condition from the first condition for each pairing, to determine whether each gene was up- or down-regulated

OCvLC.difference <- pairwise.differences[,1:2]
ALAvLC.difference <- pairwise.differences[,c(1,3)]
ALAvOC.difference <- pairwise.differences[,c(1,4)]
LAvLC.difference <- pairwise.differences[,c(1,5)]
LAvOC.difference <- pairwise.differences[,c(1,6)]
ALAvLA.difference <- pairwise.differences[,c(1,7)]
# combines the gene ID with the difference (positive or negative) for all ANOVA differentially expressed genes for each pairwise comparison

list.OCvLC <- OCvLC.difference[OCvLC.info,]
list.ALAvLC <- ALAvLC.difference[ALAvLC.info,]
list.ALAvOC <- ALAvOC.difference[ALAvOC.info,]
list.LAvLC <- LAvLC.difference[LAvLC.info,]
list.LAvOC <- LAvOC.difference[LAvOC.info,]
list.ALAvLA <- ALAvLA.difference[ALAvLA.info,]
# reduces the genes to just those determined to be differentially expressed in the pairwise conditions by t tests

OCvLC.updown <- data.frame(Upregulated=(NA) ,Downregulated=(NA))
ALAvLC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
ALAvOC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
LAvLC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
LAvOC.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))
ALAvLA.updown <- data.frame(Upregulated=(NA), Downregulated=(NA))

for (i in 1:nrow(list.OCvLC)){
	if (list.OCvLC[i,2]>0){
		OCvLC.updown[i,1] <- list.OCvLC[i,1]
	}
	if (list.OCvLC[i,2]<0){
		OCvLC.updown[i,2] <- list.OCvLC[i,1]
	}
}

for (i in 1:nrow(list.ALAvLC)){
	if (list.ALAvLC[i,2]>0){
		ALAvLC.updown[i,1] <- list.ALAvLC[i,1]
	}
	if (list.ALAvLC[i,2]<0){
		ALAvLC.updown[i,2] <- list.ALAvLC[i,1]
	}
}

for (i in 1:nrow(list.ALAvOC)){
	if (list.ALAvOC[i,2]>0){
		ALAvOC.updown[i,1] <- list.ALAvOC[i,1]
	}
	if (list.ALAvOC[i,2]<0){
		ALAvOC.updown[i,2] <- list.ALAvOC[i,1]
	}
}

for (i in 1:nrow(list.LAvLC)){
	if (list.LAvLC[i,2]>0){
		LAvLC.updown[i,1] <- list.LAvLC[i,1]
	}
	if (list.LAvLC[i,2]<0){
		LAvLC.updown[i,2] <- list.LAvLC[i,1]
	}
}

for (i in 1:nrow(list.LAvOC)){
	if (list.LAvOC[i,2]>0){
		LAvOC.updown[i,1] <- list.LAvOC[i,1]
	}
	if (list.LAvOC[i,2]<0){
		LAvOC.updown[i,2] <- list.LAvOC[i,1]
	}
}

for (i in 1:nrow(list.ALAvLA)){
	if (list.ALAvLA[i,2]>0){
		ALAvLA.updown[i,1] <- list.ALAvLA[i,1]
	}
	if (list.ALAvLA[i,2]<0){
		ALAvLA.updown[i,2] <- list.ALAvLA[i,1]
	}
}
# puts gene ID into the up- or down-regulated column appropriately for each pariwise comparison

OCvLC.up <- na.omit(OCvLC.updown[,1])
ALAvLC.up <- na.omit(ALAvLC.updown[,1])
ALAvOC.up <- na.omit(ALAvOC.updown[,1])
LAvLC.up <- na.omit(LAvLC.updown[,1])
LAvOC.up <- na.omit(LAvOC.updown[,1])
ALAvLA.up <- na.omit(ALAvLA.updown[,1])

OCvLC.down <- na.omit(OCvLC.updown[,2])
ALAvLC.down <- na.omit(ALAvLC.updown[,2])
ALAvOC.down <- na.omit(ALAvOC.updown[,2])
LAvLC.down <- na.omit(LAvLC.updown[,2])
LAvOC.down <- na.omit(LAvOC.updown[,2])
ALAvLA.down <- na.omit(ALAvLA.updown[,2])
# gets rid of NA values and separates the up- and down-regulated gene IDs into separate lists


### Write Up- and Down-Regulated Genes to Files for FunNet ###

dir.create(paste(subsubdir.DE.pair, "Up Regulated Genes", sep="/"))
subsubdir.DE.pairU <- paste(subsubdir.DE.pair, "Up Regulated Genes/", sep="/")

write.table(OCvLC.up, paste(subsubdir.DE.pairU, "OCvLC Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLC.up, paste(subsubdir.DE.pairU, "ALAvLC Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvOC.up, paste(subsubdir.DE.pairU, "ALAvOC Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvLC.up, paste(subsubdir.DE.pairU, "LAvLC Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvOC.up, paste(subsubdir.DE.pairU, "LAvOC Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLA.up, paste(subsubdir.DE.pairU, "ALAvLA Upregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)

dir.create(paste(subsubdir.DE.pair, "Down Regulated Genes", sep="/"))
subsubdir.DE.pairD <- paste(subsubdir.DE.pair, "Down Regulated Genes/", sep="/")

write.table(OCvLC.down, paste(subsubdir.DE.pairD, "OCvLC Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLC.down, paste(subsubdir.DE.pairD, "ALAvLC Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvOC.down, paste(subsubdir.DE.pairD, "ALAvOC Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvLC.down, paste(subsubdir.DE.pairD, "LAvLC Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvOC.down, paste(subsubdir.DE.pairD, "LAvOC Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLA.down, paste(subsubdir.DE.pairD, "ALAvLA Downregulated Genes.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Determine Number of DE, Up- and Down-Regulated Genes #####

total.num <- c(length(OCvLC.genes), length(ALAvLC.genes), length(ALAvOC.genes), length(LAvLC.genes), length(LAvOC.genes), length(ALAvLA.genes))
up.num <- c(length(OCvLC.up), length(ALAvLC.up), length(ALAvOC.up), length(LAvLC.up), length(LAvOC.up), length(ALAvLA.up))
down.num <- c(length(OCvLC.down), length(ALAvLC.down), length(ALAvOC.down), length(LAvLC.down), length(LAvOC.down), length(ALAvLA.down))

DE.numbers <- data.frame(Comparison=(pairwise.comparisons), Total=(total.num), UpRegulated=(up.num), DownRegulated=(down.num))

write.table(DE.numbers, paste(subsubdir.DE, "Differentially Expressed Gene Numbers.txt", sep=""), quote=F, row.names=F, sep="\t")

ANOVA.DE.number <- c("Total DE Genes Determined by ANOVA", length(p.05))

write.table(c("\n", ANOVA.DE.number), paste(subsubdir.DE, "Differentially Expressed Gene Numbers.txt", sep=""), append=T, quote=F, row.names=F, col.names=F, sep="\t")
