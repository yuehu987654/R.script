library(devtools)
library(ggplot2)
library(mvMORPH)
library(RPANDA)
library(rstan)
library(rstanarm)
rstan_options(auto_write = TRUE)
library(RColorBrewer)
source("./ABDOMEN.R")
tree <- read.tree("./199.tre")
table <- t(read.table("./199.txt", header=TRUE, sep="\t",row.names = 1))
name <- "run_bacterial_orders" 
code_path <- getwd()
detection_threshold <- 1e-05
seed <- 3
mean_prior_logY <- 0
sd_prior_logY <- 2
nb_cores <- 4
chains <-  4
warmup <-  1000
iter <-  2000
fit_summary <- ABDOMEN(tree, table, name, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)
library(ggimage)
ABDOMEN_process_output(tree, table, name, fit_summary)
original_lambda <- ABDOMEN_extract_lambda(tree, table, fit_summary) 
original_lambda
ABDOMEN_extract_Z0(tree, table, fit_summary)
R_matrices <- ABDOMEN_extract_R(tree, table, fit_summary) 
R_matrices$R 
R_matrices$R_lower_bound
R_matrices$R_upper_bound
R_matrices$R_signif
name_random <- "run_112"
seed <- 100
set.seed(seed)
table_random <- table[sample(tree$tip.label),]
rownames(table_random) <- rownames(table)

fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random, 
                              code_path = code_path,
                              detection_threshold = detection_threshold, seed = seed, 
                              mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                              nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)

ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut)
nb_permutations <- 100
list_lambda_permutations <- c()
for (seed in 1:nb_permutations){
  set.seed(seed)
  name_random <- paste0("run_Cetartiodactyla_bacterial_orders_permutation_",seed)
  table_random <- table[sample(tree$tip.label),] # randomly permutates all the Cetartiodactyla species:
  rownames(table_random) <- rownames(table)
  fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random, 
                                code_path = code_path,
                                detection_threshold = detection_threshold, seed = seed, 
                                mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                                nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)
  list_lambda_permutations <- rbind(list_lambda_permutations, ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut))
}
length(which(list_lambda_permutations[,1]>=original_lambda[1]))/nb_permutations 