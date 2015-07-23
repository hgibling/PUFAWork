load("~/Desktop/SummerWork/PUFA All.Rdata")
library(gage)
library(stringr)
library(KEGGREST)


### Generate Gene Sets for Rat KEGG Pathways

all.genes <- as.numeric(row.names(no.duplicates))

rat.gene.sets <- kegg.gsets(species="rat", id.type="entrez")
kegg.sets <- rat.gene.sets$kg.sets


# GMT File

set.names <- names(kegg.sets)
split.names <- str_split(set.names, " ", n=2)

for (i in 1:length(set.names)){
	mat <- matrix(nrow=1, c(split.names[[i]][2], split.names[[i]][1], kegg.sets[[i]]))
	write.table(mat, "Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
}


## Get Missing KEGG Pathways

# Manually compared list of Rat KEGG Pathway IDs obtained above to list on KEGG website: 25 pathways missing. 7 pathways are not able to have genes retrieved from keggGet--need to figure out solution

missing <- c("rno00220", "rno04014", "rno04015", "rno04022", "rno04024", "rno04068", "rno04071", "rno04152", "rno04261", "rno04550", "rno04611", "rno04750", "rno04919", "rno04921", "rno04922", "rno04923", "rno05230", "rno05231")

for (i in 1:length(missing)){
	pathway.info <- keggGet(missing[i])
	full.name <- pathway.info[[1]]$NAME
	name <- gsub(" -.*", "", full.name)
	full.genes <- pathway.info[[1]]$GENE
	get <- seq(1, length(full.genes), 2)
	genes <- full.genes[get]
	mat <- matrix(nrow=1, c(name, missing[i], genes))
	write.table(mat, "Rat KEGG Sets.gmt", sep="\t", quote=F, row.names=F, col.names=F, append=T)
}


# GCT File

gct.first <- "#1.2"
gct.second <- c(nrow(no.duplicates), ncol(no.duplicates))

gsea.data <- data.frame(NAME=rownames(no.duplicates), Description=rep("Gene", nrow(no.duplicates)), no.duplicates)

write.table(gct.first, "GSEA for All.gct", quote=F, col.names=F, row.names=F)
write.table(gct.second, "GSEA for All.gct", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(gsea.data, "GSEA for All.gct", quote=F, row.names=F, append=T, sep="\t")


# CLS File

cls.first <- matrix(nrow=1, (c(ncol(no.duplicates), length(pca.conditions), 1)))
cls.second <- matrix(nrow=1, c("#", pca.conditions))
cls.third <- matrix(nrow=1, c(as.vector(conditions)))

write.table(cls.first, "GSEA for All.cls", quote=F, col.names=F, row.names=F, sep="\t")
write.table(cls.second, "GSEA for All.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(cls.third, "GSEA for All.cls", quote=F, col.names=F, row.names=F, append=T, sep="\t")



