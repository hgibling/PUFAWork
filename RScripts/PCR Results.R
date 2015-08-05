library(tidyr)
library(dplyr)
library(plyr)

setwd("~/Desktop/PCR Results")

# Prep the data by removing unnecessary columns, identifying triplicates with a standard deviation of more than 0.2, removing the furthest replicate from the median, and recalculating the triplicate mean

gene.prep <- function(filename){
	file <- read.csv(filename)
	remove.excess <- file[, c(5:8)]
	remove.excess <- filter(remove.excess, !Cq=="N/A")
	remove.excess$Cq <- as.numeric(remove.excess$Cq)
	grouped <- group_by(remove.excess, Sample)
	filtered.above <- filter(grouped, Cq.Std..Dev>0.2)
	filtered.below <- filter(grouped, Cq.Std..Dev<=0.2) %>%
		select(Sample, Gene.Cq.Mean=Cq.Mean)
	remove.furthest <- dplyr::mutate(filtered.above, Median.Cq=median(Cq)) %>%
		dplyr::mutate(Distance=abs(Median.Cq-Cq)) %>%
		filter(!Distance==max(Distance)) %>%
		select(Sample, Cq) %>%
		dplyr::mutate(Gene.Cq.Mean=mean(Cq)) %>%
		select(-Cq)
	combined <- bind_rows(remove.furthest, filtered.below)
	duplicates <- which(duplicated(combined$Sample)==T)
	sample.means <- arrange(combined[-duplicates,], Sample)
	wells.removed <- nrow(file)-nrow(combined)
	print(c("Number of wells removed:", wells.removed))
	return(sample.means)
}

housekeeping.prep <- function(filename){
	output <- gene.prep(filename)
	colnames(output) <- c("Sample", "Housekeeping.Cq.Mean")
	return(output)
}

Rn18s.cDNA1 <- housekeeping.prep("new cDNA1 Rn18s.csv")
# file name must be in quotes
Adam9.cDNA1 <- gene.prep("new cDNA1 Adam9.csv")
Col3a1.cDNA1 <- gene.prep("new cDNA1 Col3a1.csv")
Col15a1.cDNA1 <- gene.prep("new cDNA1 Col15a1.csv")
Lgi1.cDNA1 <- gene.prep("new cDNA1 Lgi1.csv")
Lyz2.cDNA1 <- gene.prep("new cDNA1 Lyz2.csv")
Pdgfd.cDNA1 <- gene.prep("new cDNA1 Pdgfd.csv")
Vwa7.cDNA1 <- gene.prep("new cDNA1 Vwa7.csv")

Rn18s.cDNA2 <- housekeeping.prep("new cDNA2 Rn18s.csv")
Angptl2.cDNA2 <- gene.prep("new cDNA2 Angptl2.csv")
Angptl4.cDNA2 <- gene.prep("new cDNA2 Angptl4.csv")
Hmcn1.cDNA2 <- gene.prep("new cDNA2 Hmcn1.csv")
Igf1.cDNA2 <- gene.prep("new cDNA2 Igf1.csv")
Nrp1.cDNA2 <- gene.prep("new cDNA2 Nrp1.csv")
Ntn4.cDNA2 <- gene.prep("new cDNA2 Ntn4.csv")
St3gal2.cDNA2 <- gene.prep("new cDNA2 St3gal2.csv")
Wnt16.cDNA2 <- gene.prep("new cDNA2 Wnt16.csv")


# Normalize the Cq values of the gene of interest for each sample by subtracting the Cq values from the housekeeping gene, then find the mean Cq for each condition

normalized.gene <- function(gene, housekeeping){
	compare <- join(gene, housekeeping, by="Sample")
	normalized <- mutate(compare, Normalized.Cq=Gene.Cq.Mean-Housekeeping.Cq.Mean) %>%
		select(Sample, Normalized.Cq)
	separated <- separate(normalized, col=Sample, into=c("Sample", "Replicate"))
	summarized <- group_by(separated, Sample) %>%
		dplyr::summarize(Average.Cq=mean(Normalized.Cq))
		return(summarized)
}

Adam9.normalized <- normalized.gene(Adam9.cDNA1, Rn18s.cDNA1)
Col3a1.normalized <- normalized.gene(Col3a1.cDNA1, Rn18s.cDNA1)
Col15a1.normalized <- normalized.gene(Col15a1.cDNA1, Rn18s.cDNA1)
Lgi1.normalized <- normalized.gene(Lgi1.cDNA1, Rn18s.cDNA1)
Lyz2.normalized <- normalized.gene(Lyz2.cDNA1, Rn18s.cDNA1)
Pdgfd.normalized <- normalized.gene(Pdgfd.cDNA1, Rn18s.cDNA1)
Vwa7.normalized <- normalized.gene(Vwa7.cDNA1, Rn18s.cDNA1)

Angptl2.normalized <- normalized.gene(Angptl2.cDNA2, Rn18s.cDNA2)
Angptl4.normalized <- normalized.gene(Angptl4.cDNA2, Rn18s.cDNA2)
Hmcn1.normalized <- normalized.gene(Hmcn1.cDNA2, Rn18s.cDNA2)
Igf1.normalized <- normalized.gene(Igf1.cDNA2, Rn18s.cDNA2)
Nrp1.normalized <- normalized.gene(Nrp1.cDNA2, Rn18s.cDNA2)
Ntn4.normalized <- normalized.gene(Ntn4.cDNA2, Rn18s.cDNA2)
St3gal2.normalized <- normalized.gene(St3gal2.cDNA2, Rn18s.cDNA2)
Wnt16.normalized <- normalized.gene(Wnt16.cDNA2, Rn18s.cDNA2)


# Determine the fold change for each pairwise comparison of conditions by using the delta delta Ct method
# 2^ because every Cq difference equals a doubling of the strand, and -1 because lower Cq means there is more product to be detected earlier and thus the gene is present in greater amounts

fold.change <- function(norm.gene, genename){
	OCvLC.fc <- 2^((-1)*(norm.gene[4,2]-norm.gene[3,2]))
	ALAvLC.fc <- 2^((-1)*(norm.gene[1,2]-norm.gene[3,2]))
	ALAvOC.fc <- 2^((-1)*(norm.gene[1,2]-norm.gene[4,2]))
	LAvLC.fc <- 2^((-1)*(norm.gene[2,2]-norm.gene[3,2]))
	LAvOC.fc <- 2^((-1)*(norm.gene[2,2]-norm.gene[4,2]))
	ALAvLA.fc <- 2^((-1)*(norm.gene[1,2]-norm.gene[2,2]))
	results <- data.frame(genename, OCvLC.fc, ALAvLC.fc, ALAvOC.fc, LAvLC.fc, LAvOC.fc, ALAvLA.fc)
	colnames(results) <- c("GeneName", "OCvLC", "ALAvLC", "ALAvOC", "LAvLC", "LAvOC", "ALAvLA")
	for (i in 2:ncol(results)){
		if (results[1,i]<1){
			results[1,i] <- -1/(results[1,i])
		}
	}
	return(results)
}

Adam9.fc <- fold.change(Adam9.normalized, "Adam9")
Angptl2.fc <- fold.change(Angptl2.normalized, "Angptl2")
Angptl4.fc <- fold.change(Angptl4.normalized, "Angptl4")
Col15a1.fc <- fold.change(Col15a1.normalized, "Col15a1")
Col3a1.fc <- fold.change(Col3a1.normalized, "Col3a1")
Hmcn1.fc <- fold.change(Hmcn1.normalized, "Hmcn1")
Igf1.fc <- fold.change(Igf1.normalized, "Igf1")
Lgi1.fc <- fold.change(Lgi1.normalized, "Lgi1")
Lyz2.fc <- fold.change(Lyz2.normalized, "Lyz2")
Nrp1.fc <- fold.change(Nrp1.normalized, "Nrp1")
Ntn4.fc <- fold.change(Ntn4.normalized, "Ntn4")
Pdgfd.fc <- fold.change(Pdgfd.normalized, "Pdgfd")
St3gal2.fc <- fold.change(St3gal2.normalized, "St3gal2")
Vwa7.fc <- fold.change(Vwa7.normalized, "Vwa7")
Wnt16.fc <- fold.change(Wnt16.normalized, "Wnt16")

PCR.results <- rbind(Adam9.fc, Angptl2.fc, Angptl4.fc, Col15a1.fc, Col3a1.fc, Hmcn1.fc, Igf1.fc, Lgi1.fc, Lyz2.fc, Nrp1.fc, Ntn4.fc, Pdgfd.fc, St3gal2.fc, Vwa7.fc, Wnt16.fc)

write.table(PCR.results, "~/Desktop/PCR fold changes.csv", row.names=F, col.names=T, quote=F, sep=",")


# Test for significance of fold changes

ttest.fc <- function(gene, housekeeping, genename){
	compare <- join(gene, housekeeping, by="Sample")
	normalized.fc <- mutate(compare, Normalized.Cq=Gene.Cq.Mean-Housekeeping.Cq.Mean) %>%
		select(Sample, Normalized.Cq)
	mat <- matrix(nrow=1, c(normalized.fc[,2]))
	colnames(mat) <- normalized.fc[,1]
	rownames(mat) <- genename
	return(mat)
}

Adam9.ttest <- ttest.fc(Adam9.cDNA1, Rn18s.cDNA1, "Adam9")
Angptl2.ttest <- ttest.fc(Angptl2.cDNA2, Rn18s.cDNA2, "Angptl2")
Angptl4.ttest <- ttest.fc(Angptl4.cDNA2, Rn18s.cDNA2, "Angptl4")
Col15a1.ttest <- ttest.fc(Col15a1.cDNA1, Rn18s.cDNA1, "Col15a1")
Col3a1.ttest <- ttest.fc(Col3a1.cDNA1, Rn18s.cDNA1, "Col3a1")
Hmcn1.ttest <- ttest.fc(Hmcn1.cDNA2, Rn18s.cDNA2, "Hmcn1")
Igf1.ttest <- ttest.fc(Igf1.cDNA2, Rn18s.cDNA2, "Igf1")
Lgi1.ttest <- ttest.fc(Lgi1.cDNA1, Rn18s.cDNA1, "Lgi1")
Lyz2.ttest <- ttest.fc(Lyz2.cDNA1, Rn18s.cDNA1, "Lyz2")
Nrp1.ttest <- ttest.fc(Nrp1.cDNA2, Rn18s.cDNA2, "Nrp1")
Ntn4.ttest <- ttest.fc(Ntn4.cDNA2, Rn18s.cDNA2, "Ntn4")
Pdgfd.ttest <- ttest.fc(Pdgfd.cDNA1, Rn18s.cDNA1, "Pdgfd")
St3gal2.ttest <- ttest.fc(St3gal2.cDNA2, Rn18s.cDNA2, "St3gal2")
Vwa7.ttest <- ttest.fc(Vwa7.cDNA1, Rn18s.cDNA1, "Vwa7")
Wnt16.ttest <- ttest.fc(Wnt16.cDNA2, Rn18s.cDNA2, "Wnt16")

matbind <- rbind(Adam9.ttest, Angptl2.ttest, Angptl4.ttest, Col15a1.ttest, Col3a1.ttest, Hmcn1.ttest, Igf1.ttest, Lgi1.ttest, Nrp1.ttest, Ntn4.ttest, Pdgfd.ttest, St3gal2.ttest, Vwa7.ttest, Wnt16.ttest)

normalized.conditions <- gl(4, 8, 32, label=c("ALA", "LA", "LC", "OC"))

pcr.ttest.results <- matrix(nrow=nrow(matbind), ncol=6)
colnames(pcr.ttest.results) <- pairwise.comparisons
rownames(pcr.ttest.results) <- rownames(matbind)

for (i in 1:nrow(matbind)){
	p <- pairwise.t.test(matbind[i,], normalized.conditions, p.adj="none")
	pv <- p$p.value
	pcr.ttest.results[i,1] <- pv[3,3]
	pcr.ttest.results[i,2] <- pv[2,1]
	pcr.ttest.results[i,3] <- pv[3,1]
	pcr.ttest.results[i,4] <- pv[2,2]
	pcr.ttest.results[i,5] <- pv[3,2]
	pcr.ttest.results[i,6] <- pv[1,1]
}

pcr.results.sig <- PCR.results

for (i in 1:nrow(pcr.ttest.results)){
	for (j in 1:ncol(pcr.ttest.results)){
		if (pcr.ttest.results[i,j] > 0.05){
			pcr.results.sig[i,j+1] <- NA
		}
	}
}

write.table(pcr.results.sig, "~/Desktop/Significant PCR fold changes.csv", row.names=F, col.names=T, quote=F, sep=",")
