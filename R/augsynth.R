################################################################################
## Main functions for augmented synthetic controls Method
################################################################################

augsynth <- function(formula, unit, time, t_int, data) {

    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula)[[2]]
    trt <- terms(formula)[[3]]
    wide <- format_data(outcome, trt, unit, time, t_int, data)
    synth_data <- do.call(format_synth, wide)
    
    ## fit augsynth
    asyn <- fit_ridgeaug_formatted(wide, synth_data)
    
    ##format output
    class(asyn) <- "augsynth"
    return(asyn)
}
