#' @importFrom stats as.formula lm predict qnorm quantile rgamma rmultinom rnorm
utils::globalVariables(c("time", "val", "post", "weight", ".", "Time",
                         "Estimate", "Std.Error", "Level", "last_time",
                         "is_avg", "label", "Outcome", "unit", "obs",
                         "lambdas", "errors_se",
                         "upper_bound", "lower_bound",
                         ":=", "ATT", "ID", "Tx", "Yhat", "Yobs",
                         "ever_Tx", "id", "impact", "inf", "p_val",
                         "raw_average", "rstat", "sdo", "trt",
                         "trt_status", "trt_unit", "tx"))