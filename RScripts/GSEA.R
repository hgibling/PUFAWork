### All 4 Conditions ###

first.line <- "#1.2"
second.line <- c(nrow(no.duplicates), ncol(no.duplicates))

gsea.data <- data.frame(NAME=rownames(no.duplicates), Description=rep("Gene", nrow(no.duplicates)), no.duplicates)

write.table(first.line, "GSEA for All.gct", quote=F, col.names=F, row.names=F)
write.table(second.line, "GSEA for All.gct", quote=F, col.names=F, row.names=F, append=T, sep="\t")
write.table(gsea.data, "GSEA for All.gct", quote=F, row.names=F, append=T, sep="\t")

