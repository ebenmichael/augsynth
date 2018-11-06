# ents
[![Build Status](https://travis-ci.org/ebenmichael/augsynth.svg?branch=master)](https://travis-ci.org/ebenmichael/augsynth) [![Project Status: WIP  Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)


## Overview
This package implements the Augmented Synthetic Control Method (ASCM).

## Installation
To install this package, first ensure that `devtools` is installed with

```
install.packages("devtools")
```

then install the package from GitHub with

```
devtools::install_github("ebenmichael/augsynth")
```

## Basic usage
To get started, use a panel dataset with an `outcome` measure, a `treatment` indicator, a `unit` indicator, a `time` variable, and an intervention time `t_int`. Then run


```
asyn <- augsynth(outcome ~ trt, unit, time, t_int, data)
```
