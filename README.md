# ents
[![Build Status](https://travis-ci.org/ebenmichael/ents.svg?branch=master)](https://travis-ci.org/ebenmichael/ents) [![Project Status: WIP  Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)


## Overview
This package implements methods for maximum entropy synthetic controls using multiple outcomes.

## Installation
To install this package, first ensure that `devtools` is installed with

```
install.packages("devtools")
```

then install the package from GitHub with

```
devtools::install_github("ebenmichael/ents")
```

## Basic usage

The following requires a  data frame `outcomes` which at a minimum has the following columns:

- `outcome_id`: ID or name of the outcome
- `time`: Time that outcome was measured
- `outcome`: Value of the outcome
- `treated`: Whether the unit is treated

and a data frame `metadata` which at a minimum has the following columns:

- `t_int`: Time of intervention/treatment

Four methods of imputing a synthetic control are implemented:

- `get_synth`: The Abadie, Diamond, Hainmueller (2010) synthetic controls estimator, fit using [https://cran.r-project.org/web/packages/Synth/Synth.pdf](Synth)
- `get_l2_entropy`: The maximum entropy synthetic controls estimator
- `get_dr`: The maximum entropy synthetic controls estimator augmented with a linear outcome model
- `get_ipw`: An IPW estimator fit with regularized logistic regression


Each function takes in `outcomes` and `metadata` along with hyper-parameters, and returns as list including a dataframe with the synthetic control added and the synthetic control weights.
