library(pca3d)

raw.3d <- pca3d(pca.values.raw$x, components=1:3, col=pca.colors, show.plane=T, show.labels=T)
pre.3d <- pca3d(pca.values.preprocessed$x, components=1:3, col=pca.colors, show.plane=T, show.labels=T)

