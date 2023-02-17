## ---------------------------
##
## Script name: tidyClusters
##
## Purpose of script: 
##  This Rscript takes an JSON file containing the clusters of CoVizu as input, 
##    tidy it up for better visualization based on the width and sampling dates of each variant,
##    and outputs the tidied JSON file
##
## Author: Yiran Shao, yshao242@uwo.ca
##
## Date Created: 2022-05-18
##
## Last edited: 2022-07-12
## ---------------------------

## Use packages
library("purrr")      # For mapping
library("rjson")      # For processing json files
library("data.tree")  # For custom tree
# library("rbenchmark") # For benchmark
## ---------------------------

## Define functions
replace_null <- function(x){
  # Helper function to replace null from edgesList
  # 
  # Args:
  #   x: Input list in which we want to replace null values
  x <- map(x, ~ replace(.x, is.null(.x), NA_character_))
  map(x, 
      ~ if (is.list(.x)){
        replace_null(.x)
      }
      else .x
  )
}


get_date <- function(nodesList){
  # Helper function to get sequenced dates
  # 
  # Args:
  #   nodesList: Input nodes list in which sequenced dates are found
  # 
  # Returns:
  #   nodeDateStack: Output dataframe with sequence names in column ind and sequenced dates in column values
  nodeDates <- as.data.frame(matrix(nrow = length(nodesList), ncol = 4))
  colnames(nodeDates) <- c("nodeID", "firstDate", "lastDate", "length")
  nodeDates$nodeID <- names(nodesList)
  for (nodesNbr in 1:length(nodesList)){
    currNodeLength <- length(nodesList[[nodesNbr]])
    if (currNodeLength > 0){
      nodeDates[nodesNbr, 2] <- as.Date(nodesList[[nodesNbr]][[1]][1], "%Y-%m-%d")
      nodeDates[nodesNbr, 3] <- as.Date(nodesList[[nodesNbr]][[currNodeLength]][1], "%Y-%m-%d")
      nodeDates[nodesNbr, 4] <- nodeDates[nodesNbr, 3] - nodeDates[nodesNbr, 2]
    }
  }
  return(nodeDates)
}


rearrangeEdges <- function(clusterNbr, rawJSON){
  # Rearrange the edges based on ascending order of the sequencing dates by dfs
  #
  # Args:
  #   clusterNbr:	The index of the cluster being reordered
  #								Integer, must be > 1 and < length(rawJSON)
  #   rawJSON:		A list showing the input raw JSON file
  #								List of lists representing each clusters
  #
  # Returns:
  #   timeOrderedJSON: The reordered list showing the output JSON file
  
  ## Extract node and edges list for the current cluster
  nodesList <- rawJSON[[clusterNbr]][["nodes"]]
  edgesList <- rawJSON[[clusterNbr]][["edges"]]
  
  if (length(edgesList) != 0){
    ## Build name-date vector
    nodeDates <- get_date(nodesList)
    
    ## Convert edgesList to dataframe edgesDf
    edgesList <- replace_null(edgesList)
    edgesDf <- as.data.frame(matrix(unlist(edgesList), nrow=length(edgesList), byrow=TRUE))
    
    ## Add date info to the edges
    nodeDates <- merge(x = edgesDf, y = nodeDates, 
                       by.x = c("V2"), by.y = c("nodeID"))
    
    ## Order edges based on length and ascending dates
    nodeDates <- nodeDates[order(nodeDates$length, nodeDates$firstDate),]
    
    # Create a tree
    thisTree <- as.Node(nodeDates, mode = "network")
    # print(thisTree, "firstDate", "length")
    
    # Save traversal result
    travserResults <- data.frame(firstDate = thisTree$Get("firstDate"))
    travserResults$nodes <- rownames(travserResults)
    
    ## Reorder edges
    edgesDf <- edgesDf[match(travserResults$nodes, edgesDf$V2),]
    edgesDf <- edgesDf[rowSums(is.na(edgesDf)) != ncol(edgesDf),]
    
    orderList <- as.numeric(rownames(edgesDf))
    
    ## Save reordered edges
    timeOrderedEdgeList <- edgesList[orderList]
    
    # Edit the cluster currently being processed
    rawJSON[[clusterNbr]][["edges"]] <- timeOrderedEdgeList
  }
  return(rawJSON)
}


main <- function(){
  # Main function
  
  # Get arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # First argument as the path to cluster file
  inputFilename <- args[1]
  rawJSON <- fromJSON(file = inputFilename)
  
  # Process all clusters
  clusterVec <- 1:length(rawJSON)
  for (cluster in clusterVec) {
    rawJSON <- rearrangeEdges(cluster, rawJSON)
    gc()  # Force R to release memory it is no longer using
  }
  
  # Save the output cluster file
  outputFilename <- args[2]
  outputJSON <- toJSON(rawJSON)
  write(outputJSON, file = outputFilename)
}
## ---------------------------

## Program starts
main()
