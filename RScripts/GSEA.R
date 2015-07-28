load("~/Desktop/SummerWork/PUFA All.Rdata")
setwd("~/Desktop")
library(gage)
library(stringr)
library(KEGGREST)
library(dplyr)


### Generate Gene Sets for Rat KEGG Pathways

all.genes <- as.numeric(row.names(no.duplicates))

rat.gene.sets <- kegg.gsets(species="rat", id.type="entrez")
kegg.sets <- rat.gene.sets$kg.sets


### Write Files for GSEA Analysis

# GMT File

set.names <- names(kegg.sets)
split.names <- str_split(set.names, " ", n=2)

for (i in 1:length(set.names)){
	mat <- matrix(nrow=1, c(split.names[[i]][2], split.names[[i]][1], kegg.sets[[i]]))
	write.table(mat, "~/Desktop/SummerWork/Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
}


## Get Missing KEGG Pathways

# Manually compared list of Rat KEGG Pathway IDs obtained above to list on KEGG website: 2 pathways missing.

missing <- c("rno00220", "rno04923")

for (i in 1:length(missing)){
	pathway.info <- keggGet(missing[i])
	full.name <- pathway.info[[1]]$NAME
	name <- gsub(" -.*", "", full.name)
	full.genes <- pathway.info[[1]]$GENE
	get <- seq(1, length(full.genes), 2)
	genes <- full.genes[get]
	mat <- matrix(nrow=1, c(name, missing[i], genes))
	write.table(mat, "~/Desktop/SummerWork/Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
}


# GCT File

gct.first <- "#1.2"
gct.second <- c(nrow(no.duplicates), ncol(no.duplicates))

gsea.data <- data.frame(NAME=rownames(no.duplicates), Description=rep("Gene", nrow(no.duplicates)), no.duplicates)

write.table(gct.first, "~/Desktop/SummerWork/GSEA Expression.gct", quote=F, col.names=F, row.names=F)
write.table(gct.second, "~/Desktop/SummerWork/GSEA Expression.gct", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(gsea.data, "~/Desktop/SummerWork/GSEA Expression.gct", quote=F, row.names=F, append=T, sep="\t")


# CLS File

cls.first <- matrix(nrow=1, (c(ncol(no.duplicates), length(pca.conditions), 1)))
cls.second <- matrix(nrow=1, c("#", pca.conditions))
cls.third <- matrix(nrow=1, c(as.vector(conditions)))

write.table(cls.first, "~/Desktop/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, sep="\t")
write.table(cls.second, "~/Desktop/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(cls.third, "~/Desktop/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")


### Analyzing Results

significant.gene.sets <- function(filename, comparison){
	filename <- read.csv(paste("~/Desktop/GSEA Results/", filename, ".csv", sep=""), stringsAsFactors=F)
	file.des <- arrange(filename, FDR.q.val)
	sig.pos <- which(file.des$FDR.q.val<0.25) 
	sig.sets <- file.des$NAME[sig.pos]
	mat <- matrix(nrow=1, c(comparison, length(sig.sets), sig.sets))
	return(mat)
}

OCvLC.up.GSEA <- significant.gene.sets("OCvLC", "OCvLC upregulated")
OCvLC.down.GSEA <- significant.gene.sets("LCvOC", "OCvLC downregulated")
ALAvLC.up.GSEA <- significant.gene.sets("ALAvLC", "ALAvLC upregulated")
ALAvLC.down.GSEA <- significant.gene.sets("LCvALA", "ALAvLC downregulated")
ALAvOC.up.GSEA <- significant.gene.sets("ALAvOC", "ALAvOC upregulated")
ALAvOC.down.GSEA <- significant.gene.sets("OCvALA", "ALAvOC downregulated")
LAvLC.up.GSEA <- significant.gene.sets("LAvLC", "LAvLC upregulated")
LAvLC.down.GSEA <- significant.gene.sets("LCvLA", "LAvLC downregulated")
LAvOC.up.GSEA <- significant.gene.sets("LAvOC", "LAvOC upregulated")
LAvOC.down.GSEA <- significant.gene.sets("OCvLA", "LAvOC downregulated")
ALAvLA.up.GSEA <- significant.gene.sets("ALAvLA", "ALAvLA upregulated")
ALAvLA.down.GSEA <- significant.gene.sets("LAvALA", "ALAvLA downregulated")

table <- t(rbind.fill.matrix(OCvLC.up.GSEA, OCvLC.down.GSEA, ALAvLC.up.GSEA, ALAvLC.down.GSEA, ALAvOC.up.GSEA, ALAvOC.down.GSEA, LAvLC.up.GSEA, LAvLC.down.GSEA, LAvOC.up.GSEA, LAvOC.down.GSEA, ALAvLA.up.GSEA, ALAvLA.down.GSEA))

write.table(table, "Significant GSEA Results.txt", quote=F, row.names=F, col.names=F, sep="\t")

