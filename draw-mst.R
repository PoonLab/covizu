setwd('~/git/covizu/')

# load cluster info
clusters <- read.csv('data/clusters.info.csv', header=T)
names(clusters) <- c('label', 'parent', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})

# load first component of split MST
files <- Sys.glob('mst/*.edgelist.csv')
f <- files[1]

edges <- read.csv(f, header=T)

# assign earliest collection date of child node to edges in component
edges$coldate <- sapply(edges$child, function(x) {
  min(clusters$coldate[clusters$accession==x])
})
edges$coldate <- as.Date(edges$coldate, origin='1970-01-01')


nodes <- data.frame(accession=unique(edges$parent, edges$child))
nodes$coldate <- sapply(edges$parent, function(x) {
  min(clusters$coldate[clusters$accession==x])
})
nodes$coldate <- as.Date(nodes$coldate, origin='1970-01-01')


par(mar=c(2,0,0,0))
plot(NA, xlim=range(c(edges$c.date, edges$p.date)), 
     ylim=c(1, nrow(edges)),
     bty='n', xaxt='n'
     )
for (e in edges) {
  segments(x0=e$p.date, x)
}

