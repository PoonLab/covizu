require(igraph)
require(jsonlite)
require(Rtsne)

# open TN93 distance matrix
cat("loading TN93 distance matrix\n")
tn93 <- read.csv('data/variants.tn93.txt', skip=1, header=F)
stopifnot(nrow(tn93) == ncol(tn93))


# apply hierarchical clustering to distance matrix
cat("hierarchical clustering\n")

# direct clustering
# hc <- hclust(as.dist(tn93), method='complete')
# clusters <- cutree(hc, h=0.002)
# hist(log10(table(clusters)), breaks=20, col='grey', border='white')

# cluster t-stochastic neighbor embedding
set.seed(1)
res <- Rtsne(tn93, is_distance=TRUE, verbose=TRUE)
hc <- hclust(dist(res$Y))
clusters <- cutree(hc, h=10)

# use regex to extract headers from FASTA
headers <- rep(NA, times=nrow(tn93))
con <- file('data/variants.fa', open='r')
i <- 1
while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
  if (grepl("^>", line)) {
    headers[i] <- gsub("^>(.+)", "\\1", line)
    i <- i + 1
  }
}
close(con)
stopifnot(length(headers) == nrow(tn93))

# extract sample collection dates from headers
dates <- as.Date(sapply(headers, function(x) {
  strsplit(x, "\\|")[[1]][3]
  }))

# identify the earliest sample (Wuhan, IPBCAMS-WH-01)
root <- which.min(dates)


# open CSV with SARS-COV-2 genome variant information
variants <- read.csv('data/variants.csv')
stopifnot(all(is.element(variants$cluster, headers)))


result <- lapply(1:max(clusters), function(i) {
  cat ('.')
  
  # extract cluster indices to map to headers vector
  idx <- as.integer(which(clusters==i))
  
  if (length(idx)==1) {
    list(nodes=headers[idx], edges=NA)
  } 
  else {
    # find earliest variant in cluster that is closest to root
    min.dist <- min(tn93[root, idx])
    candidates <- headers[idx[which(tn93[root, idx] == min.dist)]]
    subroot <- candidates[which.min(sapply(candidates, function(x) { 
      as.Date(strsplit(x, "\\|")[[1]][3])
      }))]
    
    # generate minimum spanning tree
    mx <- as.matrix(tn93[idx, idx])
    colnames(mx) <- headers[idx]
    g <- graph.adjacency(mx, weighted=TRUE)
    g.mst <- igraph::mst(g, algorithm='prim')
    
    # traverse MST and export node and edge lists
    el <- get.edgelist(g.mst)
    
    traverse <- function(node, parent, edgelist, edges=c()) {
      if (!is.na(parent)) {
        edges <- c(edges, parent, node)  
      }
      
      # get local edges
      temp <- el[apply(el, 1, function(e) is.element(node, e)), ]
      temp <- unique(as.vector(temp))
      children <- temp[!is.element(temp, c(node, parent))]
      
      for (child in children) {
        edges <- traverse(child, node, edgelist, edges)
      }
      return(edges)
    }
    edges <- traverse(subroot, NA, edgelist)

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

    # shorten edge list to accession numbers only
    edges <- gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", edges)
    edges <- matrix(edges, ncol=2, byrow=TRUE)

    list(nodes=nodes, edges=edges)
  }
})
cat ('\nwriting JSON file\n')
write(toJSON(result, pretty=TRUE), file="data/clusters.json")

