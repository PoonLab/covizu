setwd('~/git/covizu/')

# load cluster info
clusters <- read.csv('data/clusters.info.csv', header=T, stringsAsFactors = FALSE)
names(clusters) <- c('label', 'node.name', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})
clusters$desc <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][1]
})


#' plot.mst
#' Generate a beadplot given an edgelist (CSV) describing the 
#' variants in a subtree of a minimum spanning tree.
#' 
#' @param f:  path to CSV summarizing an MST cluster
#' @param threshold:  integer, minimum number of instances for a given 
#'                    variant to be labelled in the plot
#' @param mar:  margin argument for par()
#' @param col1:  colour for horizontal (variant) line segments
#' @param col2:  colour for vertical (link) line segments
plot.mst <- function(f, threshold=5, mar=c(2,6,0,1)+0.1, 
                     col1='black', col2='slategray2', ...) {
  edges <- read.csv(f, header=T, stringsAsFactors = FALSE)
  
  # assign earliest collection date of child node to edges in component
  edges$coldate <- sapply(edges$child, function(x) {
    min(clusters$coldate[clusters$node.name==x])
  })
  edges$coldate <- as.Date(edges$coldate, origin='1970-01-01')
  
  nodes <- data.frame(label=unique(
    c(edges$parent, edges$child)
  ))
  
  # match node names in cluster to the sample collection date as 
  # well as the coldates of related nodes
  temp <- sapply(nodes$label, function(a) {
    c(edges$coldate[edges$parent==a],
      clusters$coldate[clusters$label==a])
  })
  nodes$mindate <- as.Date(sapply(temp, min), origin='1970-01-01')
  nodes$maxdate <- as.Date(sapply(temp, max), origin='1970-01-01')
  nodes$count <- sapply(temp, length)

  # sort nodes by pre-order traversal
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
  nodes <- nodes[match(res, nodes$label), ]
  #nodes <- nodes[order(nodes$mindate), ]
  nodes$y <- 1:nrow(nodes)
  
  # horizontal range
  x <- c(nodes$mindate, nodes$maxdate)
  
  # prepare plot region
  par(mar=mar)
  plot(NA, xlim=range(x), ylim=c(1, nrow(nodes)), 
       ylab='', bty='n', xaxt='n', yaxt='n', ...)
  axis(side=1, at=pretty(x), label=strftime(pretty(x), format='%b %d'), 
       cex.axis=0.8, mgp=c(4,0.5,0))
  
  # draw vertical segments representing edges in minimum spanning tree
  for (j in 1:nrow(edges)) {
    e <- edges[j, ]
    n1 <- nodes[nodes$label==e$parent, ]
    n2 <- nodes[nodes$label==e$child, ]
    arrows(x0=e$coldate, y0=n1$y, y1=n2$y, length=0.05, col=col2)  
  }
  
  # draw horizontal segments representing genome variants
  par(xpd=NA)
  for (i in 1:nrow(nodes)) {
    n <- nodes[i,]
    segments(x0=n$mindate, x1=n$maxdate, y0=n$y, col=col1)
    
    # draw "beads"
    cluster <- clusters[clusters$label==n$label,]
    temp <- table(cluster$coldate)
    reg <- split(cluster$region, cluster$coldate)
    
    points(x=as.Date(names(temp)), y=rep(n$y, length(temp)), 
           cex=sqrt(temp), pch=21, 
           bg=sapply(reg, function(x) {
             ifelse('Canada' %in% x, 'red', 'white')
             })
    )
    if (n$count > threshold) {
      label <- gsub("hCoV-19/([^/]+)/([^/]+)/[0-9]+", "\\1/\\2", 
                    clusters$desc[clusters$label==n$label][1])
      text(min(n$mindate)-0.5, n$y, label=label, cex=0.5, adj=1)
    }
  }
  par(xpd=FALSE)
}

# batch processing
files <- Sys.glob('mst/component-*.edgelist.csv')
for (f in files) {
  # look ahead to see number of edges
  edges <- read.csv(f)
  nodes <- c(as.character(edges$parent), as.character(edges$child))
  n.nodes <- length(unique(nodes))
  
  #pdf(file=gsub("\\.edgelist\\.csv", ".pdf", f), 
  #    width=5, height=n.nodes/10)
  png(file=gsub("\\.edgelist\\.csv", ".png", f), 
      width=5*150, height=n.nodes/10*150, res=150)
  plot.mst(f, mar=c(2,4,1,1), threshold=0)
  dev.off()
}

