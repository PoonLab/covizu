setwd('~/git/covizu/')

clusters <- read.csv('data/clusters.info.csv', header=T)
#names(clusters) <- c('label', 'parent', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})

phy <- read.tree('2020-04-15-0001_treetime/clean.nwk')
ladderize(phy, right=F)
# some clusters were dropped due to incomplete collection dates
clusters <- clusters[is.element(clusters$accession, phy$tip.label), ]


#ti <- max(pd$nodes$x) - pd$nodes$x 
#axis(side=1, at=max(pd$nodes$x) - pretty(pd$nodes$x), labels=pretty(ti))

est.dates <- read.table('2020-04-15-0001_treetime/dates.tsv', sep='\t')
names(est.dates) <- c('node.name', 'isodate', 'float.date')
est.dates$date <- as.Date(est.dates$isodate)

#all(pd$nodes$label[1:Ntip(phy)] == names(by.acc)[index])
index <- match(pd$nodes$label, est.dates$node.name)
pd$nodes$min.date <- est.dates$date[index]

# convert from ISO date to tree height (fractional years)
fit <- lm(x ~ min.date, data=pd$nodes)
date.to.x <- function(dt) {
  as.numeric(predict(fit, newdata=data.frame(min.date=dt)))
}

# associate tips with cluster info
by.acc <- split(clusters, clusters$accession)


# TODO: re-calculate vertical location, sorting by internal node time

require(ggfree)
pd <- as.phyloData(phy)

L <- tree.layout(phy)
index <- match(L$nodes$label, names(by.acc))

setwd('~/git/covizu')
png(file='vignettes/tree.png', width=4*150, height=4*150, res=150)

plot(L, label='n', mar=c(0,0,0,0), xlim=c(0, 0.6), lwd=1)
#axis(side=1)
#for (i in 1:length(by.acc)) {
#  stopifnot(names(by.acc)[index[i]] == L$nodes$label[i])
#  cluster <- by.acc[[index[i]]]
#  xr <- range(cluster$coldate)
#  if (xr[1] < xr[2]) {
#    segments(x0=date.to.x(xr[1]), x1=date.to.x(xr[2]), 
#             y0=L$nodes$y[i], lwd=2)  
#  }
#}
dev.off()

