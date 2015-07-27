load("~/Desktop/SummerWork/PUFA All.Rdata")
setwd("~/Desktop")
library(gage)
library(stringr)
library(KEGGREST)


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
	write.table(mat, "/SummerWork/Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
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
	write.table(mat, "/SummerWork/Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
}


# GCT File

gct.first <- "#1.2"
gct.second <- c(nrow(no.duplicates), ncol(no.duplicates))

gsea.data <- data.frame(NAME=rownames(no.duplicates), Description=rep("Gene", nrow(no.duplicates)), no.duplicates)

write.table(gct.first, "/SummerWork/GSEA Expression.gct", quote=F, col.names=F, row.names=F)
write.table(gct.second, "/SummerWork/GSEA Expression.gct", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(gsea.data, "/SummerWork/GSEA Expression.gct", quote=F, row.names=F, append=T, sep="\t")


# CLS File

cls.first <- matrix(nrow=1, (c(ncol(no.duplicates), length(pca.conditions), 1)))
cls.second <- matrix(nrow=1, c("#", pca.conditions))
cls.third <- matrix(nrow=1, c(as.vector(conditions)))

write.table(cls.first, "/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, sep="\t")
write.table(cls.second, "/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(cls.third, "/SummerWork/GSEA Classes.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")


### Analyzing Results




