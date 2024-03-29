% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{conformal_inf}
\alias{conformal_inf}
\title{Conformal inference procedure to compute p-values and point-wise confidence intervals}
\usage{
conformal_inf(
  ascm,
  alpha = 0.05,
  stat_func = NULL,
  type = "iid",
  q = 1,
  ns = 1000,
  grid_size = 50
)
}
\arguments{
\item{ascm}{Fitted `augsynth` object}

\item{alpha}{Confidence level}

\item{stat_func}{Function to compute test statistic}

\item{type}{Either "iid" for iid permutations or "block" for moving block permutations; default is "block"}

\item{q}{The norm for the test static `((sum(x ^ q))) ^ (1/q)`}

\item{ns}{Number of resamples for "iid" permutations}

\item{grid_size}{Number of grid points to use when inverting the hypothesis test}
}
\value{
List that contains:
        \itemize{
         \item{"att"}{Vector of ATT estimates}
         \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
         \item{"se"}{Standard error, always NA but returned for compatibility}
         \item{"lb"}{Lower bound of 1 - alpha confidence interval}
         \item{"ub"}{Upper bound of 1 - alpha confidence interval}
         \item{"p_val"}{p-value for test of no post-treatment effect}
         \item{"alpha"}{Level of confidence interval}
        }
}
\description{
Conformal inference procedure to compute p-values and point-wise confidence intervals
}
