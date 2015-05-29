# genes of interest were obtained as follows:
##		Entrez gene IDs from the ANOVA list of DE genes were run through Batch Entrez
##		the 'related database' option was chosen for RefSeq proteins
##		the list of RefSeq protein IDs was downloaded in fasta format ordered by assession number
##		the XP_ ids were removed (predicted proteins)
##		the protein fasta sequences were run through several signal peptide prediction servers
##		a running list of genes/proteins with predicted signal peptide sequences was created
##		the final list of genes of interest were those with signal peptide sequences as predicted by Signal-BLAST, TOPCONS, and SignalP (two levels of sensitivity)

library(dplyr)
load("~/Desktop/SummerWork/PUFA.Rdata")

interest <- read.csv("~/Desktop/SummerWork/GenesofInterest.csv", header=T)


# remove gene ids with protein isoforms

dupes <- which(duplicated(interest[,1]))
interest <- interest[-dupes,]
interest.tbldf <- tbl_df(interest)


# work with sensitive SignalP results

interest.sensitive.tbldf <- select(interest.tbldf, -SignalP4.1)


# genes predicted by all 3 methods to have signal peptide regions

sensitive3tbldf <- filter(interest.sensitive.tbldf, SignalBLAST==T, TOPCONS==T, SignalP.sensitive==T)
sensitive3 <- as.data.frame(sensitive3tbldf)
sensitive3ids <- sort(sensitive3[,1])


# genes predicted by at least one method to have signal peptide regions

interest.sensitive <- interest[,-4]
sensitive.ids <- sort(interest.sensitive[,1])


# determine which genes of interest are DE in pairwise comparisons

interest.in.pairwise <- function(interest, pairwise.list){
	pos <- pairwise.list %in% interest
	pairwise.list[which(pos==T)]
	#return interest.list
}

# DE predicted by all 3 methods

OCvLC.interest.s3 <- interest.in.pairwise(sensitive3ids, OCvLC.genes)
ALAvLC.interest.s3 <- interest.in.pairwise(sensitive3ids, ALAvLC.genes)
ALAvOC.interest.s3 <- interest.in.pairwise(sensitive3ids, ALAvOC.genes)
LAvLC.interest.s3 <- interest.in.pairwise(sensitive3ids, LAvLC.genes)
LAvOC.interest.s3 <- interest.in.pairwise(sensitive3ids, LAvOC.genes)
ALAvLA.interest.s3 <- interest.in.pairwise(sensitive3ids, ALAvLA.genes)


# up-regulated predicted by all 3 methods

OCvLC.interest.s3.up <- interest.in.pairwise(sensitive3ids, OCvLC.up)
ALAvLC.interest.s3.up <- interest.in.pairwise(sensitive3ids, ALAvLC.up)
ALAvOC.interest.s3.up <- interest.in.pairwise(sensitive3ids, ALAvOC.up)
LAvLC.interest.s3.up <- interest.in.pairwise(sensitive3ids, LAvLC.up)
LAvOC.interest.s3.up <- interest.in.pairwise(sensitive3ids, LAvOC.up)
ALAvLA.interest.s3.up <- interest.in.pairwise(sensitive3ids, ALAvLA.up)


# down-regulated predicted by all 3 methods

OCvLC.interest.s3.down <- interest.in.pairwise(sensitive3ids, OCvLC.down)
ALAvLC.interest.s3.down <- interest.in.pairwise(sensitive3ids, ALAvLC.down)
ALAvOC.interest.s3.down <- interest.in.pairwise(sensitive3ids, ALAvOC.down)
LAvLC.interest.s3.down <- interest.in.pairwise(sensitive3ids, LAvLC.down)
LAvOC.interest.s3.down <- interest.in.pairwise(sensitive3ids, LAvOC.down)
ALAvLA.interest.s3.down <- interest.in.pairwise(sensitive3ids, ALAvLA.down)


# DE predicted by at least one method

OCvLC.interest.s <- interest.in.pairwise(sensitive.ids, OCvLC.genes)
ALAvLC.interest.s <- interest.in.pairwise(sensitive.ids, ALAvLC.genes)
ALAvOC.interest.s <- interest.in.pairwise(sensitive.ids, ALAvOC.genes)
LAvLC.interest.s <- interest.in.pairwise(sensitive.ids, LAvLC.genes)
LAvOC.interest.s <- interest.in.pairwise(sensitive.ids, LAvOC.genes)
ALAvLA.interest.s <- interest.in.pairwise(sensitive.ids, ALAvLA.genes)


# up-regulated predicted by at least one method

OCvLC.interest.s.up <- interest.in.pairwise(sensitive.ids, OCvLC.up)
ALAvLC.interest.s.up <- interest.in.pairwise(sensitive.ids, ALAvLC.up)
ALAvOC.interest.s.up <- interest.in.pairwise(sensitive.ids, ALAvOC.up)
LAvLC.interest.s.up <- interest.in.pairwise(sensitive.ids, LAvLC.up)
LAvOC.interest.s.up <- interest.in.pairwise(sensitive.ids, LAvOC.up)
ALAvLA.interest.s.up <- interest.in.pairwise(sensitive.ids, ALAvLA.up)


# down-regulated predicted by at least one method

OCvLC.interest.s.down <- interest.in.pairwise(sensitive.ids, OCvLC.down)
ALAvLC.interest.s.down <- interest.in.pairwise(sensitive.ids, ALAvLC.down)
ALAvOC.interest.s.down <- interest.in.pairwise(sensitive.ids, ALAvOC.down)
LAvLC.interest.s.down <- interest.in.pairwise(sensitive.ids, LAvLC.down)
LAvOC.interest.s.down <- interest.in.pairwise(sensitive.ids, LAvOC.down)
ALAvLA.interest.s.down <- interest.in.pairwise(sensitive.ids, ALAvLA.down)


# write DE numbers to table

dir.create(paste(subsubdir.DE.pair, "DE Genes of Interest", sep="/"))
subsubdir.DE.interest<-paste(subsubdir.DE.pair,"DE Genes of Interest/",sep="/")

interest.total.s3 <- c(length(OCvLC.interest.s3), length(ALAvLC.interest.s3), length(ALAvOC.interest.s3), length(LAvLC.interest.s3), length(LAvOC.interest.s3), length(ALAvLA.interest.s3))
interest.up.s3 <- c(length(OCvLC.interest.s3.up), length(ALAvLC.interest.s3.up), length(ALAvOC.interest.s3.up), length(LAvLC.interest.s3.up), length(LAvOC.interest.s3.up), length(ALAvLA.interest.s3.up))
interest.down.s3 <- c(length(OCvLC.interest.s3.down), length(ALAvLC.interest.s3.down), length(ALAvOC.interest.s3.down), length(LAvLC.interest.s3.down), length(LAvOC.interest.s3.down), length(ALAvLA.interest.s3.down))

interest.numbers.s3 <- data.frame(Comparison=pairwise.comparisons, Total=interest.total.s3, UpRegulated=interest.up.s3, Downregulated=interest.down.s3)

interest.s3.tally <- c("Total DE genes of interest predicted to have a peptide signal sequence by all 3 methods", length(sensitive3ids))

interest.total.s <- c(length(OCvLC.interest.s), length(ALAvLC.interest.s), length(ALAvOC.interest.s), length(LAvLC.interest.s), length(LAvOC.interest.s), length(ALAvLA.interest.s))
interest.up.s <- c(length(OCvLC.interest.s.up), length(ALAvLC.interest.s.up), length(ALAvOC.interest.s.up), length(LAvLC.interest.s.up), length(LAvOC.interest.s.up), length(ALAvLA.interest.s.up))
interest.down.s <- c(length(OCvLC.interest.s.down), length(ALAvLC.interest.s.down), length(ALAvOC.interest.s.down), length(LAvLC.interest.s.down), length(LAvOC.interest.s.down), length(ALAvLA.interest.s.down))

interest.numbers.s <- data.frame(Comparison=pairwise.comparisons, Total=interest.total.s, UpRegulated=interest.up.s, Downregulated=interest.down.s)

interest.s.tally <- c("Total DE genes of interest predicted to have a peptide signal sequence by at least one method", length(sensitive.ids))

write.table(interest.numbers.s3, paste(subsubdir.DE.interest, "Genes of Interest Numbers 3 Methods.txt", sep=""), quote=F, row.names=F, sep="\t")
write.table(c("\n", interest.s3.tally), paste(subsubdir.DE.interest, "Genes of Interest Numbers 3 Methods.txt", sep=""), append=T, quote=F, row.names=F, col.names=F, sep="\t")

write.table(interest.numbers.s, paste(subsubdir.DE.interest, "Genes of Interest Numbers One or More.txt", sep=""), quote=F, row.names=F, sep="\t")
write.table(c("\n", interest.s.tally), paste(subsubdir.DE.interest, "Genes of Interest Numbers One or More.txt", sep=""), append=T, quote=F, row.names=F, col.names=F, sep="\t")
