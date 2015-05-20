together <- c(seq(from=1, to=32, by=4), seq(from=2, to=32, by=4), seq(from=3, to=32, by=4), seq(from=4, to=32, by=4))

tog.colors <- c(rep("red", 8), rep("green", 8), rep("cyan", 8), rep("purple", 8))

boxplot(raw, range=1.5, col=pca.colors, xlab="Array", ylab="Log Probe Intensity", main="Raw Log Probe Intensity")

hist(raw, col=pca.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=pca.colors, lty=1)
hist(raw, col=hist.colors, lty=1, xlab="Log Intensity", ylab="Density", main="Raw Density Estimation")
legend("topright", inset=0.01, cex=0.75, c(condition.names), col=hist.colors, lty=1)

black3 <- rep("black", 3)
black2 <- rep("black", 2)
hist.colors <- c(black2, "red", black3, "blue", black3, "green", black3, "orange", black3, "cyan", black3, "purple", black3, "brown", black3, "yellow", "black")


dev.off()

?RLE
