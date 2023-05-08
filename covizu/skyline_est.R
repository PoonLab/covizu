options(warn = -1)

library(ape)
library(rphyloxml)
library("optparse")
library(phytools)
library(LambdaSkyline)
library(data.table)


set.seed(123456)
# nboots <- 10

#Get key parameters controlling the simulations
list_of_options = list(
  make_option(c("-t", "--tree_file"), type = "character", 
              help = "phyloxml file containing the tree to be analyzed"),
  make_option(c("-d", "--dates_file"), type = "character", 
              help = "file containing dates associated with each tip"),
  make_option(c("-l", "--sequence_labels_file"), type = "character", 
              help = "file containing labels for all tips")
)

# Read in data
input_parser = OptionParser(option_list = list_of_options)
args = parse_args(input_parser)
tree = read.nexus(args$tree_file)
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
      tree <- bind.tip(tree, paste0(tip_name,"_", counter), edge.length = 0, 
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


# 
# #Function to run boostrapping on data, drops one branch of the phylogeny
# get_pop_sizes <- function(rand_tips, tree){
#   boot_tree <- drop.tip(tree, rand_tips)
#   alpha = betacoal.maxlik(boot_tree)
#   skyline = (skyline.multi.phylo(boot_tree, alpha$p1))
#   pop_sizes <- head(skyline$population.size, n = 5)
#   return(mean(pop_sizes, na.rm = TRUE))
# }
# 
# 
# ntips <- length(tree$tip.label)
# replicate_pop_size <- apply(matrix(sample.int(ntips, nboots*ceiling(ntips*0.05), replace = F), ncol = round(ntips*0.05)), 
#                             MARGIN = 1,
#                             get_pop_sizes, 
#                             tree = tree)
# replicate_pop_size <- replicate_pop_size[order(replicate_pop_size)]
# t0 <- mean(replicate_pop_size)
# lb <- replicate_pop_size[ceiling(0.025* nboots)]
# ub <- replicate_pop_size[round(0.975* nboots)]
# 
# 
# cat(paste(t0, lb,ub, sep = ","))
# 
# 
# 