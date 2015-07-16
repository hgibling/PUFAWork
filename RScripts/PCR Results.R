library(tidyr)
library(dplyr)

setwd("~/Desktop")

gene.rep.cq.means <- function(filename, genename){
	file <- read.csv(filename)
	just.rep.mean <- file[,c(5,7)]
	dup.pos <- which(duplicated(just.rep.mean)==T)
	just.mean <- just.rep.mean[-c(dup.pos),]
	cq.mean.name <- paste("Av.Cq", genename, sep=".")
	colnames(just.mean) <- c("Sample", cq.mean.name)
	return(just.mean)
}

rn18s <- gene.rep.cq.means("Rn18s check.csv", "Rn18s")
# file name and gene name must be in quotes
angptl4 <- gene.rep.cq.means("Angptl4 check.csv", "Angptl4")

compare <- merge(rn18s, angptl4)
compare <- separate(compare, col=Sample, into=c("Sample", "Replicate"))
compare <- compare[,-2]
compare.just <- group_by(compare, Sample) %>% summarize(MeanCqRn18s=mean(Av.Cq.Rn18s), MeanCqAngptl4=mean(Av.Cq.Angptl4))

normalized <- mutate(compare.just, Normalized=MeanCqAngptl4-MeanCqRn18s)

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


