##################################################
## Script to evaluate estimators in a parallel way
##################################################
source("evaluate_sim.R")
## get command line arguments

args <- commandArgs(TRUE)

#if(length(args)==0){
#    print("No arguments supplied.")
#    ##supply default values
#    a = 1
#    b = c(1,1,1)
#}else{
#    for(i in 1:length(args)){
#      eval(parse(text=args[[i]]))
#    }
#}

#n_units <- args[1]
#t_total <- args[2]
#t_int <- args[3]
#d <- args[4]
#lambda <- args[5]
#corrs_min <- args[6]
#corrs_max <- args[7]
#corrs_len <- args[8]
#n_sims <- args[9]

n_units <- 100
t_total <- 25
t_int <- 15
d <- 10
lambda <- 1
corrs_min <- -.9
corrs_max <- .9
corrs_len <- 10
n_sims <- 3000

## get the environment variable for number of clusters
n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
print(n_cores)

## run simulations and get results

res <- eval_factor(n_units, t_total, t_int, d, lambda,
                   seq(corrs_min, corrs_max, length.out = corrs_len),
                   n_sims,
                   n_cores = n_cores)

## read claned basque data

#bas_outcomes <- read.csv("data/clean_basque.csv")
#bas_metadata <- read.csv("data/clean_basque_meta.csv")
#bas <- list(outcomes=bas_outcomes, metadata=bas_metadata)
#res <- eval_loocv(bas, seq(corrs_min, corrs_max, length.out = corrs_len),
#                  n_sims, n_cores = n_cores)
                           

## save file
timestamp <- format(Sys.time(), "%d-%b-%Y-%H-%M-%S")

model_name <- "cluster_multi_factor"
path <- paste("results", model_name, sep="/")

mse_name <- paste(path,
                  "_mse_",
                  timestamp,
                  ".csv", sep="")

write.csv(res$mse, mse_name, row.names=FALSE)

w_diff_name <- paste(path,
                     "_w_diff_",
                     timestamp,
                     ".csv", sep="")

write.csv(res$w_diff, w_diff_name, row.names=FALSE)

metadata_name <- paste(path,
                       "_metadata_",
                       timestamp,
                       ".csv", sep="")

write.csv(res$metadata, metadata_name, row.names=FALSE)
