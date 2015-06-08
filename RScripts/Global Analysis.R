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