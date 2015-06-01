# genes of interest were obtained as follows:
##		Entrez gene IDs from the ANOVA list of DE genes were run through Batch Entrez
##		the 'related database' option was chosen for RefSeq proteins
##		the list of RefSeq protein IDs was downloaded in fasta format ordered by assession number
##		the XP_ ids were removed (predicted proteins)
##		the protein fasta sequences were run through several signal peptide prediction servers
##		a running list of genes/proteins with predicted signal peptide sequences was created
##		the final list of genes of interest were those with signal peptide sequences as predicted by Signal-BLAST, TOPCONS, and SignalP, and had evidence of secretion from the UniProt database


load("~/Desktop/SummerWork/PUFA.Rdata")
interest <- read.csv("~/Desktop/SummerWork/GenesofInterest.csv", header=F)
interest <- c(interest[,1])

# determine which genes of interest are DE in pairwise comparisons

interest.in.pairwise <- function(interest.list, pairwise.list){
	pos <- pairwise.list %in% interest.list
	pairwise.list[which(pos==T)]
}

# differentially expressed

OCvLC.interest <- interest.in.pairwise(interest, OCvLC.genes)
ALAvLC.interest <- interest.in.pairwise(interest, ALAvLC.genes)
ALAvOC.interest <- interest.in.pairwise(interest, ALAvOC.genes)
LAvLC.interest <- interest.in.pairwise(interest, LAvLC.genes)
LAvOC.interest <- interest.in.pairwise(interest, LAvOC.genes)
ALAvLA.interest <- interest.in.pairwise(interest, ALAvLA.genes)


# up-regulated

OCvLC.interest.up <- interest.in.pairwise(interest, OCvLC.up)
ALAvLC.interest.up <- interest.in.pairwise(interest, ALAvLC.up)
ALAvOC.interest.up <- interest.in.pairwise(interest, ALAvOC.up)
LAvLC.interest.up <- interest.in.pairwise(interest, LAvLC.up)
LAvOC.interest.up <- interest.in.pairwise(interest, LAvOC.up)
ALAvLA.interest.up <- interest.in.pairwise(interest, ALAvLA.up)


# down-regulated

OCvLC.interest.down <- interest.in.pairwise(interest, OCvLC.down)
ALAvLC.interest.down <- interest.in.pairwise(interest, ALAvLC.down)
ALAvOC.interest.down <- interest.in.pairwise(interest, ALAvOC.down)
LAvLC.interest.down <- interest.in.pairwise(interest, LAvLC.down)
LAvOC.interest.down <- interest.in.pairwise(interest, LAvOC.down)
ALAvLA.interest.down <- interest.in.pairwise(interest, ALAvLA.down)


# write DE numbers to table

dir.create("~/Desktop/PUFA Microarray Analysis/All 32 Arrays/Preprocessed Data/Filtered Genes Data/Differentially Expressed Gene Lists/Pairwise DE Genes/Genes of Interest")
subsubdir.DE.interest <- "~/Desktop/PUFA Microarray Analysis/All 32 Arrays/Preprocessed Data/Filtered Genes Data/Differentially Expressed Gene Lists/Pairwise DE Genes/Genes of Interest"

interest.total <- c(length(OCvLC.interest), length(ALAvLC.interest), length(ALAvOC.interest), length(LAvLC.interest), length(LAvOC.interest), length(ALAvLA.interest))
interest.up <- c(length(OCvLC.interest.up), length(ALAvLC.interest.up), length(ALAvOC.interest.up), length(LAvLC.interest.up), length(LAvOC.interest.up), length(ALAvLA.interest.up))
interest.down <- c(length(OCvLC.interest.down), length(ALAvLC.interest.down), length(ALAvOC.interest.down), length(LAvLC.interest.down), length(LAvOC.interest.down), length(ALAvLA.interest.down))

interest.numbers <- data.frame(Comparison=pairwise.comparisons, Total=interest.total, UpRegulated=interest.up, Downregulated=interest.down)

interest.tally <- c("Total DE genes of interest predicted to have a peptide signal sequence with evidence of secretion", length(interest))

write.table(interest.numbers, paste(subsubdir.DE.interest, "Genes of Interest Numbers.txt", sep="/"), quote=F, row.names=F, sep="\t")
write.table(c("\n", interest.tally), paste(subsubdir.DE.interest, "Genes of Interest Numbers.txt", sep="/"), append=T, quote=F, row.names=F, col.names=F, sep="\t")
