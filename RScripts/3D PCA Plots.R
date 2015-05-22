library(pca3d)

raw.3d <- pca3d(pca.values.raw$x, components=1:3, col=plot.colors, show.plane=T, show.labels=T)
pre.3d <- pca3d(pca.values.preprocessed$x, components=1:3, col=plot.colors, show.plane=T, show.labels=T)

ann.3d <- pca3d(pca.annotated.genes$x, components=1:3, col=plot.colors, show.plane=T, show.labels=T)
