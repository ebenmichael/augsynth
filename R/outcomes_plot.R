

#' Plot function returning the original level of the outcome variable
#' for the treated unit and its synthetic counterfactual
#'
#' @param augsynth Augsynth object to be plotted
#' @param  measure Whether to plot the synthetic counterfactual or the
#'   raw average of donor units
#' @export
augsynth_outcomes_plot <- function(augsynth, measure = c("synth", "average")) {

  trt_index <- which(augsynth$data$trt == 1)
  df <- bind_cols(augsynth$data$X, augsynth$data$y)
  synth_unit <- t(df[-trt_index, ]) %*% augsynth$weights
  average_unit <- df[-trt_index, ] %>% colMeans()
  treated_unit <- t(df[trt_index, ])

  max_y <- max(c(synth_unit, average_unit, treated_unit))
  min_y <- min(c(synth_unit, average_unit, treated_unit))

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(aes(x = augsynth$data$time, y = treated_unit, linetype = augsynth$trt_unit))

  if ('synth' %in% measure) {
    p <- p +
      ggplot2::geom_line(aes(x = augsynth$data$time, y = synth_unit, linetype = 'Synthetic counterfactual'))
  }

  if ('average' %in% measure) {
    p <- p +
      ggplot2::geom_line(aes(x = augsynth$data$time, y = average_unit, linetype = 'Donor raw average'))
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, 140)) +
    ggplot2::labs(linetype = NULL,
                  x = augsynth$time_var,
                  y = 'Outcome') +
    ggplot2::ylim(min_y, max_y) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.75, 0.88),
                   legend.key = element_rect(fill = alpha("white", 0.5)),
                   legend.background = element_rect(fill = alpha("white", 0)))

  return(p)
}

