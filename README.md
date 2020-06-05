# augsynth: Augmented Synthetic Control Method
[![Build Status](https://travis-ci.org/ebenmichael/augsynth.svg?branch=master)](https://travis-ci.org/ebenmichael/augsynth) [![Project Status: Active  The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)





## Overview
This package implements the Augmented Synthetic Control Method (ASCM).

For a more detailed description of the main functionality check out:
- [the vignette for simultaneous adoption](https://github.com/ebenmichael/augsynth/blob/master/vignettes/singlesynth-vignette.md)
- [the vignette for staggered adoption](https://github.com/ebenmichael/augsynth/blob/master/vignettes/multisynth-vignette.md)

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
