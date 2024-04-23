require(jsonlite)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("\n\nUsage: Rscript infection_prediction.R [tree filename] [labels filename] \n",
       "[summary statistics JSON] [prediction filename],\n",
}

tree.filename <- args[1]
labels.filename <- args[2]
stats.filename <- args[3]
prediction.filename <- args[4]

increasing_mods <- readRDS('infections_increasing_model_comparisons.rds')
infections_mods <- readRDS('num_infections_model_comparisons.rds')
sum_stat_dat <- fromJSON(stats.filename)
estimate_vals <- function(models, predict_dat, exp = FALSE){
    prediction_df <- data.frame(sapply(models, predict, newdata = predict_dat, type = "response"))
    if (exp) {
        prediction_df <- exp(prediction_df)
    }
    return(prediction_df)
}

set.seed(123456)

tree = read.nexus(tree.filename)
sequence_labels = read.csv(labels.filename)

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
mean_pop_size <- mean(pop_sizes, na.rm = TRUE)

sum_stat_dat$Ne <- mean_pop_size
increasing_predict_prob <- estimate_vals(increasing_mods, sum_stat_dat)

if(!is.nan(sum_stat_dat$Ne)) {
    pred_prob <- increasing_predict_prob[which(rownames(increasing_predict_prob) == "HUNePi.1"), ]
} else {
    pred_prob <- increasing_predict_prob[which(rownames(increasing_predict_prob) == "HUPi.1"), ]
}

predicted_increase <- pred_prob > 0.5

if(predicted_increase){
    infections_predictions <- 
        estimate_vals(infections_mods, sum_stat_dat, exp = T)
    if(!is.nan(sum_stat_dat$Ne)) {
        predicted_infections <- infections_predictions[which(rownames(infections_predictions) == "HUNePi.1"), ]
    } else {
        predicted_infections <- infections_predictions[which(rownames(infections_predictions) == "HUPi.1"), ]
    }
} else {
    predicted_infections <- -1
}

con <- file(prediction.filename)
writeLines(paste(predicted_infections), con)