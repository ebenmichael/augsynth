% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multisynth_class.R
\name{predict.multisynth}
\alias{predict.multisynth}
\title{Get prediction of average outcome under control or ATT}
\usage{
\method{predict}{multisynth}(object, att = F, att_weight = NULL, bs_weight = NULL, ...)
}
\arguments{
\item{object}{Fit multisynth object}

\item{att}{If TRUE, return the ATT, if FALSE, return imputed counterfactual}

\item{att_weight}{Weights to place on individual units/cohorts when averaging}

\item{bs_weight}{Weight to perturb units by for weighted bootstrap}

\item{...}{Optional arguments}
}
\value{
Matrix of predicted post-treatment control outcomes for each treated unit
}
\description{
Get prediction of average outcome under control or ATT
}
