#############################################################
## Methods for choosing the balance in SC
#############################################################


bin_search_ <- function(eps, feasfunc) {
    #' Perform binary search for the smallest balance which is feasible
    #' @param eps Sorted list of balance tolerances to try
    #' @param feasfunc Function which returns True if feasible
    #'
    #' @return smallest tolerance which is feasible

    if(length(eps) == 1) {
        if(feasfunc(eps[1])) {
            return(eps[1])
        } else {
            return(-1)
        }
    } else {
        ## check feasibility of middle point (rounded down)
        mid <- floor(length(eps) / 2)
        isfeas <- feasfunc(eps[mid])
        if(isfeas) {
            ## if feasible, try lower half
            return(bin_search_(eps[1:mid], feasfunc))
        } else {
            ## if not feasible try upper half
            return(bin_search_(eps[(mid+1):length(eps)], feasfunc))
        }
    }
}


bin_search <- function(start, end, by, feasfunc) {
    #' Wrapper for bin_search_ which generates tolerances to search over
    #' @param start Starting value of tolerances
    #' @param end Ending value of tolerances
    #' @param by Step size of tolerances
    #' @param feasfunc Function which returns True if feasible

    eps <- seq(start, end, by)
    return(bin_search_(eps, feasfunc))
}


lexical <- function(outcomes, metadata, grp_order, outcome_col, trt_unit=1, by=.1, maxep=1) {
    #' Finds the lowest feasible tolerance in each group, holding
    #' the previous groups fixed
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param grp_order Lexical ordering of groups
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param by Step size for tolerances to dry, default: 0.1
    #' @param maxep Maximum tolerance to consider
    #'
    #' @return List with lowest tolerances lexically
    #' @export

    ## just format data once
    data_out <- format_data(outcomes, metadata, trt_unit, outcome_col)

    ## get the magnitdues for the controls in each group
    syn_data <- data_out$synth_data
    x <- syn_data$Z0
    groups <- data_out$groups
    mags <- lapply(groups, function(g) sqrt(sum((x[g,] - mean(x[g,]))^2)))

    ## reorder
    mags <- mags[grp_order]

    ## set original epsilon to infinity
    epslist <- lapply(grp_order, function(g) 10^20 * mags[[g]])
    names(epslist) <- grp_order

    ## iterate over groups, finding the lowest feasible tolerance
    for(g in grp_order) {

        ## create the feasibility function by keeping other tolerances fixed
        feasfunc <- function(ep) {
            epslist[[g]] <- ep
            feas <- suppressMessages(fit_entropy_formatted(data_out, epslist)$feasible)
            return(feas)
        }
        
        ## find the best epsilon
        minep <- bin_search(0, maxep, by, feasfunc)

        ## if it failed, then stop everything
        if(minep < 0) {
            break
        }
        ## keep that tolerance
        epslist[[g]] <- minep
    }

    return(epslist)
}


lexical_time <- function(outcomes, metadata, trt_unit=1, by=.1) {
    #' Finds the lowest feasible tolerance in each time period from
    #' most recent to oldest, holding the previous times fixed
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param by Step size for tolerances to dry, default: 0.1 * magnitude
    #'
    #' @return List with lowest tolerances lexically
    #' @export

    ## get the times
    t_int <- metadata$t_int
    times <- rev(paste((outcomes %>% filter(time < t_int) %>%
                    distinct(time) %>% select(time))$time))

    ## add a second time column to selection
    timeout <- outcomes %>%
        mutate(time2=ifelse(time < (t_int-1), paste(time), paste(t_int-1)))

    ## get the lowest feasible tolerances
    time_eps <- lexical(timeout, metadata, times, "time2", trt_unit, by)

    ## fit the control
    lex <- get_entropy(timeout, metadata, trt_unit, time_eps, "time2")

    ## get rid of the second time column
    lex$outcomes <- lex$outcomes %>% select(-time2)

    return(lex)
}



recent_group <- function(outcomes, metadata, t_past, trt_unit=1, by=.1, maxep=1) {
    #' Lexically minimizes the imbalance in two groups, recent and old
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param t_past Number of recent time periods to balance especially
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param by Step size for tolerances to dry, default: 0.1 * magnitude
    #' @param maxep Maximum tolerance to consider
    #'
    #' @return Fitted SC with lexically minimized imbalance in time groups
    #' @export

    ## add a column indicating before or after "most recent period"
    timeout <- outcomes %>%
        mutate(recent=ifelse(time < t_past, "Old", "Recent"))

    ## get the lowest feasible tolerances
    time_eps <- lexical(timeout, metadata, c("Recent", "Old"), "recent", trt_unit, by, maxep)
    
    ## fit the control
    lex <- get_entropy(timeout, metadata, trt_unit, time_eps, "recent")

    ## get rid of the second time column
    lex$outcomes <- lex$outcomes %>% select(-recent)

    return(lex)
}


sep_lasso_ <- function(outcomes, metadata, trt_unit, by) {
    #' Internal function that does the work of sep_lasso
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression)
    #' @param by Step size for tolerances to try
    #'
    #' @return List with lowest tolerances in units of standard deviation
    
    ## just format data once
    data_out <- format_data(outcomes, metadata, trt_unit)

    ## get the standard deviations of controls for each group
    syn_data <- data_out$synth_data
    x <- syn_data$Z0
    sds <- apply(x, 1, sd)
    ## set original epsilon to infinity
    epslist <- sapply(sds, function(sd) 10^20 * sd)
    
    ## create the feasibility function by changing the tolerance in units of sd
    feasfunc <- function(ep) {
        epslist <- sapply(sds, function(sd) ep * sd)
        feas <- suppressMessages(fit_entropy_formatted(data_out, epslist, lasso=TRUE)$feasible)
        return(feas)
    }
    
    ## find the best epsilon
    minep <- bin_search(0, 4, by, feasfunc)

    ## if it failed, then stop everything
    if(minep < 0) {
        stop("Failed to find a synthetic control with balance better than 4 std deviations")
    }
    ## keep that tolerance
    epslist <- lapply(sds, function(sd) minep * sd)
    return(epslist)
    
    }


sep_lasso <- function(outcomes, metadata, trt_unit=1, by=.1) {
    #' Finds the lowest feasible tolerance in units of standard deviation for
    #' all time periods and covariates
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param by Step size for tolerances to dry, default: 0.1 * magnitude
    #'
    #' @return max ent SC fit with lowest imbalance tolerance
    #' @export

    ## get the lowest global imbalance
    epslist <- sep_lasso_(outcomes, metadata, trt_unit, by)
    ## fit the SC
    lasso_sc <- get_entropy(outcomes, metadata, trt_unit, eps=epslist, lasso=TRUE)

    return(lasso_sc)
}
