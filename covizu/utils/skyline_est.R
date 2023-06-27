options(warn = -1)

library(LambdaSkyline)

set.seed(123456)

#Get key parameters controlling the simulations
list_of_options = list(
  optparse::make_option(c("-t", "--tree_file"), type = "character", 
              help = "nexus file containing the tree to be analyzed"),
  optparse::make_option(c("-l", "--sequence_labels_file"), type = "character", 
              help = "file containing labels for all tips")
)

# Read in data
input_parser = optparse::OptionParser(option_list = list_of_options)
args = optparse::parse_args(input_parser)
tree = ape::read.nexus(args$tree_file)
sequence_labels = read.csv(args$sequence_labels_file)
colnames(sequence_labels) = c("index", "value")

#Adjust tree to include branches of length 0 on identical sequences
tip_count = table(sequence_labels$index)
add_tip_count = data.frame(tip_count - 1)

for (tip_place_in_table in 1:nrow(add_tip_count)){
  tip_name = add_tip_count[tip_place_in_table,1]
  freq = add_tip_count[tip_place_in_table,2]
  if(freq != 0){
    for (counter in 1:freq){
      tree <- phytools::bind.tip(tree, paste0(tip_name,"_", counter), edge.length = 0, 
                           where=which(tree$tip.label == tip_name))
    }
  }
}


#Run skyline estimation

alpha = betacoal.maxlik(tree)
skyline = (skyline.multi.phylo(tree, alpha$p1))

#Output skyline estimation
pop_sizes <- head(skyline$population.size, n = 5)
cat(mean(pop_sizes, na.rm = TRUE))


