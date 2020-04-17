setwd('~/git/covizu/')

# load cluster info
clusters <- read.csv('data/clusters.info.csv', header=T)
names(clusters) <- c('label', 'node.name', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})
clusters$desc <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][1]
})

# load first component of split MST
files <- Sys.glob('mst/component-*.edgelist.csv')
f <- files[3]


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

x <- c(nodes$mindate, nodes$maxdate)

par(mar=c(2,0,0,6)+0.1)
plot(NA, xlim=range(x), ylim=c(1, nrow(nodes)), 
     bty='n', xaxt='n', yaxt='n')
axis(side=1, at=pretty(x), label=strftime(pretty(x), format='%b %d'), 
     cex.axis=0.8, mgp=c(4,0.5,0))

for (j in 1:nrow(edges)) {
  e <- edges[j, ]
  n1 <- nodes[as.character(nodes$accession)==as.character(e$parent), ]
  n2 <- nodes[as.character(nodes$accession)==as.character(e$child), ]
  arrows(x0=e$coldate, y0=n1$y, y1=n2$y, length=0.05, col='grey')  
}

for (i in 1:nrow(nodes)) {
  n <- nodes[i,]
  segments(x0=n$mindate, x1=n$maxdate, y0=n$y, col='cadetblue')
  
  # draw bubbles
  cluster <- clusters[clusters$accession==n$accession,]
  temp <- table(cluster$coldate)
  reg <- split(cluster$region, cluster$coldate)
  
  points(x=as.Date(names(temp)), y=rep(n$y, length(temp)), 
         cex=sqrt(temp), pch=21, 
         bg=sapply(reg, function(x) ifelse('Canada' %in% x, 'red', 'white'))
         )
  if (n$count > 5) {
    text(max(nodes$maxdate)+0.5, n$y, 
         label=clusters$desc[clusters$accession==n$accession][1], 
         cex=0.5, adj=0)
  }
}

