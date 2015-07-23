load("~/Desktop/SummerWork/PUFA All.Rdata")

### All 4 Conditions ###

## Generate Gene Sets for Rat KEGG Pathways
library(gage)

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



