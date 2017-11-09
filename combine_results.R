###############################################
## Script to combine outputs of simulation runs
###############################################
library(tidyverse)
library(feather)

combine_results <- function(res_list) {
    #' Combine results from different simulation runs
    #' @param res_list List of outputs from eval_sim
    #'
    #' @return Combined data as a list

    mse <- combine_dataframes(lapply(res_list, function(x) x$mse))
    w_diff <- combine_dataframes(lapply(res_list, function(x) x$w_diff))
    metadata <- combine_dataframes(lapply(res_list, function(x) x$metadata))

    return(list(mse=mse,
                metadata=metadata,
                w_diff=w_diff))
}


save_results <- function(results, model_name, directory) {
    #' Save results into a feather format
    #' @param results List of results from simulation
    #' @param model_name Name of the model simulated from
    #' @param directory Directory to save to

    df_names <- names(results)

    ## iterate over files and save
    for(i in 1:length(results)) {
        df_name <- df_names[i]
        df <- results[[i]]
        f_name <- paste(model_name, df_name, "all", sep="_")
        f_name <- paste(directory, f_name, sep="/")
        f_name <- paste(f_name, ".feather", sep="")
        ## save file
        write_feather(df, f_name)
    }
}


read_results <- function(model_name, directory) {
    #' Save results into a feather format
    #' @param model_name Name of the model simulated from
    #'
    #' @return list of results

    results = list()

    ## read all files in the directory
    files <- list.files(directory, full.names=TRUE)
    ## iterate over files and read
    i <- 1
    for(f in files) {
        if(grepl(model_name, f) & grepl("all", f)) {
            if(grepl("mse", f)) {
                df_name <- "mse"
            } else if(grepl("metadata", f)) {
                df_name <- "metadata"
            } else if (grepl("w_diff", f)) {
                df_name <- "w_diff"
            } else {
                stop(paste("Unrecognized file", f))
            }
            ## read file
            df <- as.data.frame(read_feather(f))
            results[[i]] <- df
            names(results)[i] <- df_name
            i <- i + 1
        }

    }
    return(results)
}

combine_dataframes <- function(dfs) {
    #' Row bind data frames, making sure that simulation numbers are correct
    #' @param dfs List of dataframes

    ## order by simulation number
    df_order <- order(sapply(dfs, function(df) max(df$sim_num)),
                      decreasing=TRUE)


    new_dfs <- list(dfs[[df_order[1]]])
    ## iterate over data frames and add to the simulation number
    curr_total <- max(dfs[[df_order[1]]]$sim_num)
    i <- 2
    for(idx in df_order[-1]) {
        df <- dfs[[idx]]
        n_sim <- max(df$sim_num)
        ## add to the sim number
        df$sim_num <- df$sim_num + curr_total
        ## add to the total
        curr_total <- curr_total + n_sim
        ## add to list
        new_dfs[[i]] <- df
        i <- i + 1
    }
    return(bind_rows(new_dfs))
}
