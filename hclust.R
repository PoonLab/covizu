require(igraph)
require(jsonlite)

tn93 <- read.csv('data/variants.tn93.txt', skip=1, header=F)
variants <- read.csv('data/variants.csv')

# read headers from FASTA
headers <- rep(NA, times=nrow(tn93))
con <- file('data/variants.fa', open='r')
i <- 1
while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
  if (grepl("^>", line)) {
    headers[i] <- gsub("^>(.+)", "\\1", line)
    i <- i + 1
  }
}

hc <- hclust(as.dist(tn93), method='ward.D')
clusters <- cutree(hc, h=0.002)
#hist(log10(table(clusters)), breaks=20, col='grey', border='white')

dates <- as.Date(sapply(headers, function(x) {
  strsplit(x, "\\|")[[1]][3]
  }))

# identify the earliest sample (Wuhan, IPBCAMS-WH-01)
root <- which.min(dates)

result <- lapply(1:max(clusters), function(i) {
  # find node in cluster that is closest to root
  idx <- which(clusters==i)
  if (length(idx)==1) {
    list(nodes=headers[idx], edges=NA)
  } 
  else {
    subroot <- headers[idx[which.min(tn93[root, idx])]]
    
    # generate minimum spanning tree
    mx <- as.matrix(tn93[idx, idx])
    colnames(mx) <- headers[idx]
    g <- graph.adjacency(mx, weighted=TRUE)
    g.mst <- igraph::mst(g, algorithm='prim')
    
    # traverse MST and export node and edge lists
    el <- get.edgelist(g.mst)
    
    traverse <- function(node, parent, edgelist, edges=c()) {
      # get local edges
      temp <- el[apply(el, 1, function(e) is.element(node, e)), ]
      temp <- unique(as.vector(temp))
      children <- temp[!is.element(temp, c(node, parent))]
      
      for (child in children) {
        edges <- c(edges, node, child)  
        edges <- traverse(child, node, edgelist, edges)
      }
      return(edges)
    }
    edges <- traverse(subroot, NA, edgelist)
    edges <- matrix(edges, ncol=2, byrow=TRUE)

    # store variant data
    nodes <- list()
    for (node in unique(edges)) {
      accn <- strsplit(node, "\\|")[[1]][2]
      temp <- variants[variants$cluster==node, ]
      temp$label1 <- sapply(as.character(temp$label), function(x) {
        strsplit(x, "\\|")[[1]][1]
      })
      nodes[[accn]] <- temp[c('label1', 'region', 'country', 'coldate')]
    }

    list(nodes=nodes, edges=edges)
  }
})

write(toJSON(result, pretty=TRUE), file="data/clusters.json")

# record subroots
