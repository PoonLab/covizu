setwd('~/git/covizu/')

# load cluster info
clusters <- read.csv('data/clusters.info.csv', header=T)
names(clusters) <- c('label', 'parent', 'coldate', 'region', 'country')

clusters$coldate <- as.Date(clusters$coldate)
clusters$accession <- sapply(clusters$label, function(x) {
  strsplit(as.character(x), "\\|")[[1]][2]
})

component <- read.csv('mst/cluster-0.edgelist.csv', header=T)

# assign earliest collection date of child node to edges in component
component$c.date <- sapply(component$child, function(x) {
  min(clusters$coldate[clusters$accession==x])
})
component$c.date <- as.Date(component$c.date, origin='1970-01-01')

component$p.date <- sapply(component$parent, function(x) {
  min(clusters$coldate[clusters$accession==x])
})
component$p.date <- as.Date(component$p.date, origin='1970-01-01')

component <- component[order(component$p.date, component$parent, 
                             component$c.date), ]
