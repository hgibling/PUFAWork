load("~/Desktop/SummerWork/PUFA All.Rdata")

### All 4 Conditions ###

## Convert Rat Gene IDs to Human Gene IDs (Entrez)
library(AnnotationDbi)

all.genes <- as.numeric(row.names(no.duplicates))

human.ids <- idConverter(all.genes, srcSpecies="RATNO", destSpecies="HOMSA", srcIDType="EG", destIDType="EG")

human.ids <- as.vector(human.ids)


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



