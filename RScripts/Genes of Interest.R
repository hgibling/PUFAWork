# genes of interest were obtained as follows:
##		Entrez gene IDs from the ANOVA list of DE genes were run through Batch Entrez
##		the 'related database' option was chosen for RefSeq proteins
##		the list of RefSeq protein IDs was downloaded in fasta format ordered by assession number
##		the XP_ ids were removed (predicted proteins)
##		the protein fasta sequences were run through several signal peptide prediction servers
##		a running list of genes/proteins with predicted signal peptide sequences was created
##		the final list of genes of interest were those with signal peptide sequences as predicted by Signal-BLAST, TOPCONS, and SignalP, and had evidence of secretion from the UniProt database


setwd("~/Desktop/SummerWork")
load("PUFA.Rdata")
load("PUFADE.Rdata")

interest <- read.csv("~/Desktop/SummerWork/GenesofInterest.csv", header=T)
interest.ev <- c(na.omit(interest[,1]))

# determine which genes of interest with evidence of secretion are DE in pairwise comparisons

interest.in.pairwise <- function(interest.list, pairwise.list){
	pos <- pairwise.list %in% interest.list
	pairwise.list[which(pos==T)]
}

# differentially expressed

OCvLC.interest.ev <- interest.in.pairwise(interest.ev, OCvLC.genes)
ALAvLC.interest.ev <- interest.in.pairwise(interest.ev, ALAvLC.genes)
ALAvOC.interest.ev <- interest.in.pairwise(interest.ev, ALAvOC.genes)
LAvLC.interest.ev <- interest.in.pairwise(interest.ev, LAvLC.genes)
LAvOC.interest.ev <- interest.in.pairwise(interest.ev, LAvOC.genes)
ALAvLA.interest.ev <- interest.in.pairwise(interest.ev, ALAvLA.genes)


# up-regulated

OCvLC.interest.ev.up <- interest.in.pairwise(interest.ev, OCvLC.up)
ALAvLC.interest.ev.up <- interest.in.pairwise(interest.ev, ALAvLC.up)
ALAvOC.interest.ev.up <- interest.in.pairwise(interest.ev, ALAvOC.up)
LAvLC.interest.ev.up <- interest.in.pairwise(interest.ev, LAvLC.up)
LAvOC.interest.ev.up <- interest.in.pairwise(interest.ev, LAvOC.up)
ALAvLA.interest.ev.up <- interest.in.pairwise(interest.ev, ALAvLA.up)


# down-regulated

OCvLC.interest.ev.down <- interest.in.pairwise(interest.ev, OCvLC.down)
ALAvLC.interest.ev.down <- interest.in.pairwise(interest.ev, ALAvLC.down)
ALAvOC.interest.ev.down <- interest.in.pairwise(interest.ev, ALAvOC.down)
LAvLC.interest.ev.down <- interest.in.pairwise(interest.ev, LAvLC.down)
LAvOC.interest.ev.down <- interest.in.pairwise(interest.ev, LAvOC.down)
ALAvLA.interest.ev.down <- interest.in.pairwise(interest.ev, ALAvLA.down)


# write DE numbers to table

dir.create("~/Desktop/PUFA Microarray Analysis/All 32 Arrays/Preprocessed Data/Filtered Genes Data/Differentially Expressed Gene Lists/Pairwise DE Genes/Genes of Interest")
subsubdir.DE.interest <- "~/Desktop/PUFA Microarray Analysis/All 32 Arrays/Preprocessed Data/Filtered Genes Data/Differentially Expressed Gene Lists/Pairwise DE Genes/Genes of Interest"

interest.ev.total <- c(length(OCvLC.interest.ev), length(ALAvLC.interest.ev), length(ALAvOC.interest.ev), length(LAvLC.interest.ev), length(LAvOC.interest.ev), length(ALAvLA.interest.ev))
interest.ev.up <- c(length(OCvLC.interest.ev.up), length(ALAvLC.interest.ev.up), length(ALAvOC.interest.ev.up), length(LAvLC.interest.ev.up), length(LAvOC.interest.ev.up), length(ALAvLA.interest.ev.up))
interest.ev.down <- c(length(OCvLC.interest.ev.down), length(ALAvLC.interest.ev.down), length(ALAvOC.interest.ev.down), length(LAvLC.interest.ev.down), length(LAvOC.interest.ev.down), length(ALAvLA.interest.ev.down))

interest.ev.numbers <- data.frame(Comparison=pairwise.comparisons, Total=interest.ev.total, UpRegulated=interest.ev.up, Downregulated=interest.ev.down)

interest.ev.tally <- c("Total DE genes of interest predicted to have a peptide signal sequence with evidence of secretion", length(interest.ev))

write.table(interest.ev.numbers, paste(subsubdir.DE.interest, "Genes of Interest Evidence Numbers.txt", sep="/"), quote=F, row.names=F, sep="\t")
write.table(c("\n", interest.ev.tally), paste(subsubdir.DE.interest, "Genes of Interest Evidence Numbers.txt", sep="/"), append=T, quote=F, row.names=F, col.names=F, sep="\t")



## determine which genes of interest with evidence of secretion or predicted are DE in pairwise comparisons

interest.pred <- c(na.omit(interest[,2]))


# differentially expressed

OCvLC.interest.pred <- interest.in.pairwise(interest.pred, OCvLC.genes)
ALAvLC.interest.pred <- interest.in.pairwise(interest.pred, ALAvLC.genes)
ALAvOC.interest.pred <- interest.in.pairwise(interest.pred, ALAvOC.genes)
LAvLC.interest.pred <- interest.in.pairwise(interest.pred, LAvLC.genes)
LAvOC.interest.pred <- interest.in.pairwise(interest.pred, LAvOC.genes)
ALAvLA.interest.pred <- interest.in.pairwise(interest.pred, ALAvLA.genes)


# up-regulated

OCvLC.interest.pred.up <- interest.in.pairwise(interest.pred, OCvLC.up)
ALAvLC.interest.pred.up <- interest.in.pairwise(interest.pred, ALAvLC.up)
ALAvOC.interest.pred.up <- interest.in.pairwise(interest.pred, ALAvOC.up)
LAvLC.interest.pred.up <- interest.in.pairwise(interest.pred, LAvLC.up)
LAvOC.interest.pred.up <- interest.in.pairwise(interest.pred, LAvOC.up)
ALAvLA.interest.pred.up <- interest.in.pairwise(interest.pred, ALAvLA.up)


# down-regulated

OCvLC.interest.pred.down <- interest.in.pairwise(interest.pred, OCvLC.down)
ALAvLC.interest.pred.down <- interest.in.pairwise(interest.pred, ALAvLC.down)
ALAvOC.interest.pred.down <- interest.in.pairwise(interest.pred, ALAvOC.down)
LAvLC.interest.pred.down <- interest.in.pairwise(interest.pred, LAvLC.down)
LAvOC.interest.pred.down <- interest.in.pairwise(interest.pred, LAvOC.down)
ALAvLA.interest.pred.down <- interest.in.pairwise(interest.pred, ALAvLA.down)


# write DE numbers to table

interest.pred.total <- c(length(OCvLC.interest.pred), length(ALAvLC.interest.pred), length(ALAvOC.interest.pred), length(LAvLC.interest.pred), length(LAvOC.interest.pred), length(ALAvLA.interest.pred))
interest.pred.up <- c(length(OCvLC.interest.pred.up), length(ALAvLC.interest.pred.up), length(ALAvOC.interest.pred.up), length(LAvLC.interest.pred.up), length(LAvOC.interest.pred.up), length(ALAvLA.interest.pred.up))
interest.pred.down <- c(length(OCvLC.interest.pred.down), length(ALAvLC.interest.pred.down), length(ALAvOC.interest.pred.down), length(LAvLC.interest.pred.down), length(LAvOC.interest.pred.down), length(ALAvLA.interest.pred.down))

interest.pred.numbers <- data.frame(Comparison=pairwise.comparisons, Total=interest.pred.total, UpRegulated=interest.pred.up, Downregulated=interest.pred.down)

interest.pred.tally <- c("Total DE genes of interest predicted to have a peptide signal sequence", length(interest.pred))

write.table(interest.pred.numbers, paste(subsubdir.DE.interest, "Genes of Interest Predicted Numbers.txt", sep="/"), quote=F, row.names=F, sep="\t")
write.table(c("\n", interest.pred.tally), paste(subsubdir.DE.interest, "Genes of Interest Predicted Numbers.txt", sep="/"), append=T, quote=F, row.names=F, col.names=F, sep="\t")



##### View Expression Values and Fold Changes for Each Gene for Each Condition #####

interest.ev.values.pos <- as.numeric(rownames(average.condition.values)) %in% interest.ev
interest.ev.values <- average.condition.values[which(interest.ev.values.pos==T),]

interest.pred.values.pos <- as.numeric(rownames(average.condition.values)) %in% interest.pred
interest.pred.values <- average.condition.values[which(interest.pred.values.pos==T),]

pairwise.values.ev <- data.frame(OCvLC=(interest.ev.values[,2]-interest.ev.values[,1]), ALAvLC=(interest.ev.values[,3]-interest.ev.values[,1]), ALAvOC=(interest.ev.values[,3]-interest.ev.values[,2]), LAvLC=(interest.ev.values[,4]-interest.ev.values[,1]), LAvOC=(interest.ev.values[,4]-interest.ev.values[,2]), ALAvLA=(interest.ev.values[,3]-interest.ev.values[,4]))

pairwise.ev.neg.OCvLC <- which(pairwise.values.ev[,1]<0)
pairwise.ev.neg.ALAvLC <- which(pairwise.values.ev[,2]<0)
pairwise.ev.neg.ALAvOC <- which(pairwise.values.ev[,3]<0)
pairwise.ev.neg.LAvLC <- which(pairwise.values.ev[,4]<0)
pairwise.ev.neg.LAvOC <- which(pairwise.values.ev[,5]<0)
pairwise.ev.neg.ALAvLA <- which(pairwise.values.ev[,6]<0)

pairwise.values.ev.abs <- 2^abs(pairwise.values.ev)

make.neg <- function(dataframe, column, neg.positions){
	for (i in 1:nrow(dataframe)){
		if (i %in% neg.positions){
			dataframe[i,column] <- dataframe[i,column]*(-1)
		}
	} 
	return(dataframe[,column])
}

ev.neg.OCvLC <- make.neg(pairwise.values.ev.abs, 1, pairwise.ev.neg.OCvLC)
ev.neg.ALAvLC <- make.neg(pairwise.values.ev.abs, 2, pairwise.ev.neg.ALAvLC)
ev.neg.ALAvOC <- make.neg(pairwise.values.ev.abs, 3, pairwise.ev.neg.ALAvOC)
ev.neg.LAvLC <- make.neg(pairwise.values.ev.abs, 4, pairwise.ev.neg.LAvLC)
ev.neg.LAvOC <- make.neg(pairwise.values.ev.abs, 5, pairwise.ev.neg.LAvOC)
ev.neg.ALAvLA <- make.neg(pairwise.values.ev.abs, 6, pairwise.ev.neg.ALAvLA)

pairwise.ev.fc <- data.frame(Evidence=rownames(interest.ev.values), OCvLC=ev.neg.OCvLC, ALAvLC=ev.neg.ALAvLC, ALAvOC=ev.neg.ALAvOC, LAvLC=ev.neg.LAvLC, LAvOC=ev.neg.LAvOC, ALAvLA=ev.neg.ALAvLA)

# predicted

pairwise.values.pred <- data.frame(OCvLC=(interest.pred.values[,2]-interest.pred.values[,1]), ALAvLC=(interest.pred.values[,3]-interest.pred.values[,1]), ALAvOC=(interest.pred.values[,3]-interest.pred.values[,2]), LAvLC=(interest.pred.values[,4]-interest.pred.values[,1]), LAvOC=(interest.pred.values[,4]-interest.pred.values[,2]), ALAvLA=(interest.pred.values[,3]-interest.pred.values[,4]))

pairwise.pred.neg.OCvLC <- which(pairwise.values.pred[,1]<0)
pairwise.pred.neg.ALAvLC <- which(pairwise.values.pred[,2]<0)
pairwise.pred.neg.ALAvOC <- which(pairwise.values.pred[,3]<0)
pairwise.pred.neg.LAvLC <- which(pairwise.values.pred[,4]<0)
pairwise.pred.neg.LAvOC <- which(pairwise.values.pred[,5]<0)
pairwise.pred.neg.ALAvLA <- which(pairwise.values.pred[,6]<0)

pairwise.values.pred.abs <- 2^abs(pairwise.values.pred)

pred.neg.OCvLC <- make.neg(pairwise.values.pred.abs, 1, pairwise.pred.neg.OCvLC)
pred.neg.ALAvLC <- make.neg(pairwise.values.pred.abs, 2, pairwise.pred.neg.ALAvLC)
pred.neg.ALAvOC <- make.neg(pairwise.values.pred.abs, 3, pairwise.pred.neg.ALAvOC)
pred.neg.LAvLC <- make.neg(pairwise.values.pred.abs, 4, pairwise.pred.neg.LAvLC)
pred.neg.LAvOC <- make.neg(pairwise.values.pred.abs, 5, pairwise.pred.neg.LAvOC)
pred.neg.ALAvLA <- make.neg(pairwise.values.pred.abs, 6, pairwise.pred.neg.ALAvLA)

pairwise.pred.fc <- data.frame(Evidence=rownames(interest.pred.values), OCvLC=pred.neg.OCvLC, ALAvLC=pred.neg.ALAvLC, ALAvOC=pred.neg.ALAvOC, LAvLC=pred.neg.LAvLC, LAvOC=pred.neg.LAvOC, ALAvLA=pred.neg.ALAvLA)


# write to table

write.table(pairwise.ev.fc, paste(subsubdir.DE.interest, "Genes of Interest Values.txt", sep="/"), quote=F, row.names=F, sep="\t")
write.table(pairwise.pred.fc, paste(subsubdir.DE.interest, "Genes of Interest Values.txt", sep="/"), quote=F, row.names=F, sep="\t", append=T)
# warning ok
