library(tidyr)
library(dplyr)

setwd("~/Desktop/PCR Results")

gene.prep <- function(filename){
	file <- read.csv(filename)
	remove.excess <- file[, c(5:8)]
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
# file name and gene name must be in quotes
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







OCvLC.dcq <- normalized[4,4]-normalized[3,4]
OCvLC.fc <- 2^((-1)*OCvLC.dcq)
# 2^ because every Cq difference equals a doubling of the strand, and -1 because lower Cq means there is more product to be detected earlier and thus the gene is present in greater amounts

ALAvLC.dcq <- normalized[1,4]-normalized[3,4]
ALAvLC.fc <- 2^((-1)*ALAvLC.dcq)

LAvLC.dcq <- normalized[2,4]-normalized[3,4]
LAvLC.fc <- 2^((-1)*LAvLC.dcq)


### other way ###

rn18s.summarized <- separate(rn18s, col=Sample, into=c("Sample", "Replicate"))
rn18s.summarized <- rn18s.summarized[,-2]
rn18s.just <- group_by(rn18s.summarized, Sample) %>% summarize(MeanCqRn18s=mean(Av.Cq.Rn18s))

angptl4.summarized <- separate(angptl4, col=Sample, into=c("Sample", "Replicate"))
angptl4.summarized <- angptl4.summarized[,-2]
angptl4.just <- group_by(angptl4.summarized, Sample) %>% summarize(MeanCqangptl4=mean(Av.Cq.Angptl4))

normalized2 <- angptl4.just[,2]-rn18s.just[,2]

OCvLC.dcq2 <- normalized2[4,]-normalized2[3,]
OCvLC.fc2 <- 2^((-1)*OCvLC.dcq2)

ALAvLC.dcq2 <- normalized2[1,]-normalized2[3,]
ALAvLC.fc2 <- 2^((-1)*ALAvLC.dcq2)

LAvLC.dcq2 <- normalized2[2,]-normalized2[3,]
LAvLC.fc2 <- 2^((-1)*LAvLC.dcq2)


### other other way ###

each.norm <- mutate(compare, Normalized=Av.Cq.Angptl4-Av.Cq.Rn18s)
grouped <- group_by(each.norm, Sample)
each.dcq <- summarize(grouped, MeanCqAngptl4=mean(Normalized))


