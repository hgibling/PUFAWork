### Heat Maps ###

diff.matrix <- as.matrix(diff.expressed.genes)

pdf(file="~/Desktop/Heatmap.pdf", width=20, height=20)
heatmap(diff.matrix)
dev.off()


interest.all <- c(interest.ev, interest.pred)

interest.all.values.pos <- as.numeric(rownames(diff.expressed.matrix)) %in% interest.all
interest.all.values <- diff.expressed.matrix[which(interest.all.values.pos==T),]
interest.matrix <- as.matrix(interest.all.values)

pdf(file="~/Desktop/Interest Heatmap.pdf")
heatmap(interest.matrix)
dev.off()


### Cytoscape ###

setwd("~/Desktop")

reflist1 <- read.table("Reflist1.txt", header=T)
reflist2 <- read.table("Reflist2.txt", header=T)
reflist3 <- read.table("Reflist3.txt", header=T)
reflist4 <- read.table("Reflist4.txt", header=T)
reflist5 <- read.table("Reflist5.txt", header=T)
reflist6 <- read.table("Reflist6.txt", header=T)
reflist7 <- read.table("Reflist7.txt", header=T)

reflist <- unique(rbind(reflist1, reflist2, reflist3, reflist4, reflist5, reflist6, reflist7))

expression.names.ids <- merge(reflist, remove.duplicates[,-1], by.x="GeneID", by.y="Gene")
expression.names <- expression.names.ids[,-1]
expression.ids <- expression.names.ids[,-2]

write.table(expression.names, "Expression with Gene Names.txt", row.names=F, col.names=T, sep="\t")
write.table(expression.ids, "Expression with Gene IDs.txt", row.names=F, col.names=T, sep="\t")
