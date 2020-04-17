setwd('~/git/covizu/')

# load cluster info
clusters <- read.csv('data/clusters.info.csv', header=T)
names(clusters) <- c('label', 'parent', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})

# load first component of split MST
files <- Sys.glob('mst/component-*.edgelist.csv')
f <- files[1]


edges <- read.csv(f, header=T)

# assign earliest collection date of child node to edges in component
edges$coldate <- sapply(edges$child, function(x) {
  min(clusters$coldate[clusters$accession==x])
})
edges$coldate <- as.Date(edges$coldate, origin='1970-01-01')


nodes <- data.frame(accession=unique(
  c(as.character(edges$parent), as.character(edges$child))
  ))


temp <- sapply(nodes$accession, function(a) {
  c(edges$coldate[as.character(edges$parent)==a],
    clusters$coldate[clusters$accession==a])
})
nodes$mindate <- as.Date(sapply(temp, min), origin='1970-01-01')
nodes$maxdate <- as.Date(sapply(temp, max), origin='1970-01-01')
nodes$count <- sapply(temp, length)


#nodes <- nodes[order(nodes$mindate), ]
# sort nodes in pre-order traversal
temp <- unique(edges$parent)
root <- as.character(temp[which(!is.element(temp, edges$child))])

traverse <- function(node, edges, res=c()) {
  res <- c(res, node)
  children <- edges$child[edges$parent==node]
  for (child in children) {
    res <- traverse(child, edges, res)
  }
  return(res)
}
res <- traverse(root, edges)
nodes <- nodes[match(res, nodes$accession), ]
nodes$y <- 1:nrow(nodes)


par(mar=c(2,0,0,2))
plot(NA, xlim=c(min(nodes$mindate), max(nodes$maxdate)), 
     ylim=c(1, nrow(nodes)), bty='n', xaxt='n')
for (i in 1:nrow(nodes)) {
  n <- nodes[i,]
  segments(x0=n$mindate, x1=n$maxdate, y0=n$y)
  
  # draw bubbles
  
  
  text(max(nodes$maxdate), n$y, label=n$accession, cex=0.5, adj=0)
}
for (j in 1:nrow(edges)) {
  e <- edges[j, ]
  n1 <- nodes[as.character(nodes$accession)==as.character(e$parent), ]
  n2 <- nodes[as.character(nodes$accession)==as.character(e$child), ]
  arrows(x0=e$coldate, y0=n1$y, y1=n2$y, length=0.05, col='grey')  
}

