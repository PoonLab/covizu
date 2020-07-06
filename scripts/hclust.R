require(igraph)
require(jsonlite)
require(Rtsne)

# open TN93 distance matrix
cat("loading TN93 distance matrix\n")
tn93 <- read.csv('data/variants.tn93.txt', skip=1, header=F)
tn93 <- as.matrix(tn93)
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
accns <- gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", headers)

# extract sample collection dates from headers
dates <- as.Date(sapply(headers, function(x) {
  strsplit(x, "\\|")[[1]][3]
  }))

# identify the earliest sample (Wuhan, IPBCAMS-WH-01)
root <- which.min(dates)

# compute mean pairwise distance between members of each cluster
# compute mean distance between members and the root
tab <- table(clusters)
mean.pdist <- c()
mean.rdist <- c()
for (k in 1:length(tab)) {
  clust <- as.integer(which(clusters == k))
  
  pdists <- as.matrix(tn93[clust, clust])
  mdist <- mean(pdists)
  mean.pdist <- c(mean.pdist, mdist)
  
  rdists <- as.matrix(tn93[root, clust])
  mrdist <- mean(rdists)
  mean.rdist <- c(mean.rdist, mrdist)
}

# open CSV with SARS-COV-2 genome variant information
variants <- read.csv('data/variants.csv')
stopifnot(all(is.element(variants$cluster, headers)))

#' @param node: str, label of current node variant
#' @param parent: str, label of current node's parental variant
#' @param el: str, edge list from minimum spanning tree
#' @return linearized vector of parent->child pairs
traverse <- function(node, parent, el, edges=c()) {
  if (!is.na(parent)) {
    edges <- c(edges, parent, node)  
  }
  # get adjacent to current node
  temp <- el[apply(el, 1, function(e) is.element(node, e)), ]
  temp <- unique(as.vector(temp))
  children <- temp[!is.element(temp, c(node, parent))]

  # TODO: sort children vector by genetic distance
  row <- tn93[which(headers==node), ]
  adj.dists <- row[match(children, headers)]
  adj.dates <- as.Date(gsub(".+\\|([0-9]+-[0-9]+-[0-9]+)$", "\\1", children))
  
  children <- children[order(adj.dists, adj.dates)]  # increasing
  for (child in children) {
    edges <- traverse(child, node, el, edges)
  }
  return(edges)
}

result <- list()
for (i in 1:max(clusters)) {
#result <- lapply(1:max(clusters), function(i) {
  #cat (paste0(i, '\n'))
  cat('.')
  
  # extract cluster indices to map to headers vector
  idx <- as.integer(which(clusters==i))

  # cluster mean pairwise distance
  pdist <- mean.pdist[[i]]
  
  # cluster mean root distance
  rdist <- mean.rdist[[i]]
  
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
    edges <- traverse(subroot, NA, el)

    # store variant data
    nodes <- list()
    for (node in unique(edges)) {
      accn <- strsplit(node, "\\|")[[1]][2]
      temp <- variants[variants$cluster==node, ]
      temp$label1 <- sapply(as.character(temp$label), function(x) {
        strsplit(x, "\\|")[[1]][1]
      })
      temp$accession <- sapply(as.character(temp$label), function(x) {
        strsplit(x, "\\|")[[1]][2]
      })
      nodes[[accn]] <- temp[c('label1', 'accession', 'country', 'coldate')]
    }

    # shorten edge list to accession numbers only
    edges <- gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", edges)
    edges <- matrix(edges, ncol=2, byrow=TRUE)
    
    dists <- apply(edges, 1, function(e) {
      tn93[which(accns==e[1]), which(accns==e[2])]
    })
    edges <- cbind(edges, round(dists*29903, 2))

    result[[length(result)+1]] <- list(pdist=pdist, rdist=rdist,nodes=nodes, edges=edges)
  }
}#)
cat ('\nwriting JSON file\n')
write(toJSON(result, pretty=TRUE), file="data/clusters.json")

