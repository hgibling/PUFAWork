##################################################
##################################################

# If after performing Quality Assessment it is decided all arrays should be analyzed: 

# Have the annotation file (Annotations for Rat Gene 21st.csv) on your desktop

##################################################
##################################################

library(oligo)


##### Filter Genes to Remove Those with No Annotation #####

dir.create(paste(subdir.all.preproc, "Filtered Genes Data", sep="/"))
subdir.all.preproc.filt <- paste(subdir.all.preproc, "Filtered Genes Data/", sep="/")

annotation.file <- read.csv("~/Desktop/Annotations for Rat Gene 21st.csv", header=TRUE)
annotation.matrix <- as.matrix(annotation.file)
annotations <- annotation.matrix[order(annotation.file[,1], annotation.file[,2]),]
# creates a matrix where the genes are ordered by their transcript cluster

annotated.values <- cbind(annotations, preprocessed.expression.matrix)
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
# removes repetitive columns with gene ids


### Create Reference List of All Annotated Genes for FunNet ###

reference.list <- as.numeric(rownames(no.duplicates))

write.table(reference.list, paste(subdir.all.preproc.filt, "Reference Gene List.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Images for Filtered Genes Data #####

dir.create(paste(subdir.all.preproc.filt, "Images", sep="/"))
subsubdir.images <- paste(subdir.all.preproc.filt, "Images/", sep="/")


### Boxplots ###

pdf(file=paste(subsubdir.images, "Preprossed Boxplot with Annotated Genes.pdf", sep=""), width=20, height=7)
boxplot(no.duplicates, range=1.5, col=plot.colors, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity with Annotated Genes")
dev.off()

pdf(file=paste(subsubdir.images, "Preprossed Boxplot with Annotated Genes Grouped.pdf", sep=""), width=20, height=7)
boxplot(no.duplicates[,group], range=1.5, col=group.colors, xlab="Array", ylab="Log Probe Intensity", main="Preprocessed Log Probe Intensity with Annotated Genes")
dev.off()


### PCA Plot of Filtered Genes ###

transposed.annotated.gene.values <- t(no.duplicates)
pca.annotated.genes <- prcomp(transposed.annotated.gene.values)

pdf(file=paste(subsubdir.images, "Preprocessed PCA Plot with Annotated Genes", sep=""))
par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=T)

plot(pca.annotated.genes$x, col=plot.colors, pch=20, main="Preprocessed PCA Plot with Annotated Genes")
text(pca.annotated.genes$x, pos=3, offset=0.2, labels=pca.numbers, cex=0.5)
legend("topright", inset=c(-0.15,0), c(pca.conditions), cex=0.75, col=pca.legend.colors, pch=20)
par(normal)
dev.off()

pca.summary <- summary(pca.annotated.genes)
proportion.variance <- pca.summary$importance[2:3,1:5]

pdf(file=paste(subsubdir.images, "Preprocessed Proportion of Variance PCA with Annotated Genes.pdf", sep=""))
barplot(proportion.variance, beside=T, col=c("black","gray"), main="Preprocessed Proportion of Variance of Principal Components \n with Annotated Genes", xlab="Principal Components", ylab="Proportion")
legend("topleft",inset=0.01, cex=0.75, c("Proportion of Variance", "Cumulative Proportion"), pch=15, col=c("black","gray"))
box()
dev.off()


### Hierarchical Clustering Dendogram ###

annotated.transposed <- t(no.duplicates)
annotated.distance <- dist(annotated.transposed)
annotated.sample.clusters <- hclust(annotated.distance)

pdf(file=paste(subsubdir.images, "Preprocessed Hierarchical Cluster Dendogram with Annotated Genes.pdf", sep=""))
plot(annotated.sample.clusters, main="Preprocessed Hierarchical Cluster Dendogram \n with Annotated Genes", xlab="Samples")
dev.off()


##### Rearrange Samples so Replicates are Together #####

no.duplicates <- no.duplicates[group]


##### Use ANOVA to Determine Differentially Expressed Genes #####

#code modified from Pavlidis, P. 2003 Methods 31(4):282-289

dir.create(paste(subdir.all.preproc.filt, "Differentially Expressed Gene Lists", sep="/"))
subsubdir.DE <- paste(subdir.all.preproc.filt, "Differentially Expressed Gene Lists/", sep="/")

conditions <- gl(4, 8, 32, label=c("LC", "OC", "ALA", "LA"))

anova.function <- function(x){
	data <- data.frame(conditions, x)
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

p.05 <- which(adjusted.pvalues[,1]<=0.05)

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
	ttest.results[i,1] <- pv[1,1]
	ttest.results[i,2] <- pv[2,1]
	ttest.results[i,3] <- pv[2,2]
	ttest.results[i,4] <- pv[3,1]
	ttest.results[i,5] <- pv[3,2]
	ttest.results[i,6] <- pv[3,3]
}

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
write.table(ALAvOC.genes, paste(subsubdir.DE.pairA, "ALAvOC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvLC.genes, paste(subsubdir.DE.pairA, "LAvLC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(LAvOC.genes, paste(subsubdir.DE.pairA, "LAvOC DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)
write.table(ALAvLA.genes, paste(subsubdir.DE.pairA, "ALAvLA DE Genes.txt", sep=""), quote=F, row.names=F, col.names=F)


##### Determine Up- and Down-Regulated Genes #####

average.condition.values <- data.frame(LC=rowMeans(diff.expressed.matrix[,1:8]), OC=rowMeans(diff.expressed.matrix[,9:16]), ALA=rowMeans(diff.expressed.matrix[,17:24]), LA=rowMeans(diff.expressed.matrix[,25:32]))

gene.names <- rownames(diff.expressed.genes)

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

up.down <- function(list){
	output <- data.frame(Upregulated=NA, Downregulated=NA)
	for (i in 1:nrow(list)){
		if (list[i,2]>0){
			output[i,1] <- list[i,1]
		}
		if (list[i,2]<0){
			output[i,2] <- list[i,1]
		}
	}
	return(output)
}
# puts gene ID into the up- or down-regulated column appropriately for each pariwise comparison

OCvLC.updown <- up.down(list.OCvLC)
ALAvLC.updown <- up.down(list.ALAvLC)
ALAvOC.updown <- up.down(list.ALAvOC)
LAvLC.updown <- up.down(list.LAvLC)
LAvOC.updown <- up.down(list.LAvOC)
ALAvLA.updown <- up.down(list.ALAvLA)

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
