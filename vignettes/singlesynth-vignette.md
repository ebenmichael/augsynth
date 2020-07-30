---
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Single Outcome AugSynth Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



# `augsynth`: The Augmented Synthetic Control Method


## Installation

You can install `augsynth` from github using `devtools`.


```r
## Install devtools if noy already installed
install.packages("devtools", repos='http://cran.us.r-project.org')
## Install augsynth from github
devtools::install_github("ebenmichael/augsynth")
```

## Example: Effects of the 2012 Kansas Tax Cuts 

### The data
To show the usage and features of `augsynth`, we'll use data on the impact of personal income tax cuts in Kansas that comes with the `AugSynth` package. Our interest is in estimating the effect of income tax cuts on gross state product (GSP) per capita.


```r
library(magrittr)
library(dplyr)
library(augsynth)
data(kansas)
```

The `kansas` dataset contains the GSP per capita (the outcome measure) `lngdpcapita` for all 50 states from the first quarter of 1990 to the first quarter of 2016.

To run `augsynth`, we need to include a treatment status column that indicates which region was treated and at what time. The table in `kansas` contains the column `treated` to denote this. In the original study, the second quarter of 2012 was the implementation of the tax cut in Kansas.


```r
kansas %>% 
  select(year, qtr, year_qtr, state, treated, gdp, lngdpcapita) %>% 
  filter(state == "Kansas" & year_qtr >= 2012 & year_qtr < 2013) 
#> # A tibble: 4 x 7
#>    year   qtr year_qtr state  treated    gdp lngdpcapita
#>   <dbl> <dbl>    <dbl> <chr>    <dbl>  <dbl>       <dbl>
#> 1  2012     1    2012  Kansas       0 143844        10.8
#> 2  2012     2    2012. Kansas       1 141518        10.8
#> 3  2012     3    2012. Kansas       1 138890        10.8
#> 4  2012     4    2013. Kansas       1 139603        10.8
```


### Synth
Now to find a synthetic control using the entire series of pre-intervention outcomes (and no auxiliary covariates), we can use `augsynth`. To do so we just need to give `augsynth` a formula like `outcome ~ treatment`, tell it what the unit and time variables are, optionally provide when intervention took place (the code will automatically determine this if `t_int` is not provided), and specify that we don't want to fit an outcome model


```r
library(augsynth)
syn <- augsynth(lngdpcapita ~ treated, fips, year_qtr, kansas,
                progfunc = "None", scm = T)
```

We can then look at the ATT estimates for each post-intervention time period and overall. 
We'll also see the quality of the synthetic control fit measured by the L2 distance between Kansas and its synthetic control, and the percent improvement over uniform weights.
By default, we'll also see pointwise confidence intervals using a [conformal inference procedure](https://arxiv.org/abs/1712.09089).


```r
summary(syn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "None", scm = ..2)
#> 
#> Average ATT Estimate (p Value for Joint Null):  -0.029   ( 0.331 )
#> L2 Imbalance: 0.083
#> Percent improvement from uniform weights: 79.5%
#> 
#> Avg Estimated Bias: 0.000
#> 
#> Inference type: Conformal inference
#> 
#>     Time    Estimate 95% CI Lower Bound 95% CI Upper Bound    p Value
#>  2012.25 -0.01807275        -0.04461946       0.0059457035 0.11111111
#>  2012.50 -0.04105313        -0.07012810      -0.0145064170 0.02222222
#>  2012.75 -0.03310023        -0.06217519      -0.0065535146 0.04444444
#>  2013.00 -0.01940747        -0.04595418       0.0046109870 0.11111111
#>  2013.25 -0.02915453        -0.05317298      -0.0051360786 0.04444444
#>  2013.50 -0.04626635        -0.07281306      -0.0222478978 0.02222222
#>  2013.75 -0.03181209        -0.05583054      -0.0103218972 0.02222222
#>  2014.00 -0.04476664        -0.07384161      -0.0182199278 0.02222222
#>  2014.25 -0.04267269        -0.07427592      -0.0135977204 0.02222222
#>  2014.50 -0.02940448        -0.06100771      -0.0003295131 0.04444444
#>  2014.75 -0.01849616        -0.05262764       0.0105788120 0.14444444
#>  2015.00 -0.02930564        -0.06596539       0.0048258429 0.07777778
#>  2015.25 -0.01908377        -0.05068700       0.0099911959 0.12222222
#>  2015.50 -0.02162342        -0.05575491       0.0074515445 0.11111111
#>  2015.75 -0.01856676        -0.05522650       0.0130364708 0.18888889
#>  2016.00 -0.02816962        -0.06735762       0.0084901215 0.10000000
```

It's easier to see this information visually. Below we plot the difference between Kansas and it's synthetic control. Before the tax cuts (to the left of the dashed line) we expect these to be close, and after the tax cuts we measure the effect (with point-wise confidence intervals).

<img src="figure/fig_syn-1.png" title="plot of chunk fig_syn" alt="plot of chunk fig_syn" style="display: block; margin: auto;" />

We can also compute point-wise confidence intervals using the [Jackknife+ procedure](https://arxiv.org/abs/1905.02928) by changing the `inf_type` argument, although this requires additional assumptions.

<img src="figure/fig_syn_plus-1.png" title="plot of chunk fig_syn_plus" alt="plot of chunk fig_syn_plus" style="display: block; margin: auto;" />


### Augmenting synth with an outcome model
In this example the pre-intervention synthetic control fit has an L2 imbalance of 0.083, about 20% of the imbalance between Kansas and the average of the other states. We can reduce this by _augmenting_ synth with ridge regression. To do this we change `progfunc` to `"Ridge"`. We can also choose the ridge hyper-parameter by setting `lambda`, while not specifying `lambda` will determine one through cross validation:

```r
asyn <- augsynth(lngdpcapita ~ treated, fips, year_qtr, kansas,
                progfunc = "Ridge", scm = T)
```

We can plot the cross-validation MSE when dropping pre-treatment time periods by setting `cv = T` in the `plot` function:

<img src="figure/fig_asyn_cv-1.png" title="plot of chunk fig_asyn_cv" alt="plot of chunk fig_asyn_cv" style="display: block; margin: auto;" />

By default, the CV procedure chooses the maximal value of `lambda` with MSE within one standard deviation of the minimal MSE. To instead choose the `lambda` that minizes the cross validation MSE, set `min_1se = FALSE`.


We can look at the summary and plot the results. Now in the summary output we see an estimate of the overall bias of synth; we measure this with the average amount that augmentation changes the synth estimate. Notice that the estimates become somewhat larger in magnitude, and the standard errors are tighter.

```r
summary(asyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "Ridge", scm = ..2)
#> 
#> Average ATT Estimate (p Value for Joint Null):  -0.040   ( 0.074 )
#> L2 Imbalance: 0.062
#> Percent improvement from uniform weights: 84.7%
#> 
#> Avg Estimated Bias: 0.000
#> 
#> Inference type: Conformal inference
#> 
#>     Time    Estimate 95% CI Lower Bound 95% CI Upper Bound    p Value
#>  2012.25 -0.02233432        -0.04432874        0.003043846 0.05555556
#>  2012.50 -0.04703842        -0.07580034       -0.018276493 0.02222222
#>  2012.75 -0.04252002        -0.07128194       -0.010374336 0.02222222
#>  2013.00 -0.02967474        -0.05505291       -0.004296572 0.03333333
#>  2013.25 -0.04119699        -0.06657516       -0.012435060 0.02222222
#>  2013.50 -0.05922698        -0.08798891       -0.030465057 0.02222222
#>  2013.75 -0.04466913        -0.07343106       -0.019290964 0.02222222
#>  2014.00 -0.05803598        -0.09018166       -0.025890300 0.02222222
#>  2014.25 -0.05516412        -0.09069356       -0.019634680 0.02222222
#>  2014.50 -0.04148561        -0.08039881       -0.005956177 0.03333333
#>  2014.75 -0.02923219        -0.06814538        0.006297248 0.05555556
#>  2015.00 -0.04004239        -0.08233934        0.000000000 0.05555556
#>  2015.25 -0.03045080        -0.06598023        0.001694886 0.05555556
#>  2015.50 -0.03269061        -0.07160380        0.002838827 0.05555556
#>  2015.75 -0.02895347        -0.07125042        0.009959725 0.05555556
#>  2016.00 -0.03829085        -0.08735531        0.004006102 0.05555556
```

<img src="figure/fig_asyn-1.png" title="plot of chunk fig_asyn" alt="plot of chunk fig_asyn" style="display: block; margin: auto;" />

There are also several auxiliary covariates. We can include these in the augmentation by fitting an outcome model using the auxiliary covariates. To do this we simply add the covariates into the formula after `|`. By default this will create time invariant covariates by averaging the auxiliary covariates over the pre-intervention period, dropping `NA` values. Then the lagged outcomes and the auxiliary covariates are jointly balanced by SCM and the ridge outcome model includes both.


```r
covsyn <- augsynth(lngdpcapita ~ treated | lngdpcapita + log(revstatecapita) +
                                           log(revlocalcapita) + log(avgwklywagecapita) +
                                           estabscapita + emplvlcapita,
                   fips, year_qtr, kansas,
                   progfunc = "ridge", scm = T)
```

Again we can look at the summary and plot the results.

```r
summary(covsyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "ridge", scm = ..2)
#> 
#> Average ATT Estimate (p Value for Joint Null):  -0.061   ( 0.14 )
#> L2 Imbalance: 0.054
#> Percent improvement from uniform weights: 86.6%
#> 
#> Covariate L2 Imbalance: 0.005
#> Percent improvement from uniform weights: 97.7%
#> 
#> Avg Estimated Bias: 0.000
#> 
#> Inference type: Conformal inference
#> 
#>     Time    Estimate 95% CI Lower Bound 95% CI Upper Bound    p Value
#>  2012.25 -0.02122181        -0.04426538        0.001821765 0.06666667
#>  2012.50 -0.04738718        -0.07555155       -0.014102025 0.03333333
#>  2012.75 -0.05010828        -0.08339344       -0.006581536 0.03333333
#>  2013.00 -0.04536410        -0.07352846       -0.012078937 0.03333333
#>  2013.25 -0.05507918        -0.08836434       -0.021794026 0.02222222
#>  2013.50 -0.07128399        -0.10456915       -0.032878035 0.02222222
#>  2013.75 -0.05791368        -0.09119884       -0.024628522 0.02222222
#>  2014.00 -0.08098791        -0.11939386       -0.037461166 0.02222222
#>  2014.25 -0.07794667        -0.12147342       -0.034419927 0.02222222
#>  2014.50 -0.06491981        -0.11356735       -0.021393066 0.03333333
#>  2014.75 -0.05655804        -0.11032637       -0.007910500 0.04444444
#>  2015.00 -0.07529028        -0.12393781       -0.021521943 0.03333333
#>  2015.25 -0.06264739        -0.10617414       -0.013999852 0.03333333
#>  2015.50 -0.06715143        -0.10555738       -0.018503886 0.02222222
#>  2015.75 -0.06292093        -0.10132689       -0.009152601 0.02222222
#>  2016.00 -0.07820946        -0.12173621       -0.019320335 0.02222222
```

<img src="figure/fig_covsyn-1.png" title="plot of chunk fig_covsyn" alt="plot of chunk fig_covsyn" style="display: block; margin: auto;" />

Now we can additionally fit ridge ASCM on the residuals, look at the summary, and plot the results.

```r

covsyn_resid <- augsynth(lngdpcapita ~ treated | lngdpcapita + log(revstatecapita) +
                                           log(revlocalcapita) + log(avgwklywagecapita) +
                                           estabscapita + emplvlcapita,
                   fips, year_qtr, kansas,
                   progfunc = "ridge", scm = T, lambda = asyn$lambda,
                   residualize = T)
```


```r
summary(covsyn_resid)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "ridge", scm = ..2, 
#>     lambda = ..3, residualize = ..4)
#> 
#> Average ATT Estimate (p Value for Joint Null):  -0.055   ( 0.257 )
#> L2 Imbalance: 0.067
#> Percent improvement from uniform weights: 83.4%
#> 
#> Covariate L2 Imbalance: 0.000
#> Percent improvement from uniform weights: 100%
#> 
#> Avg Estimated Bias: 0.000
#> 
#> Inference type: Conformal inference
#> 
#>     Time    Estimate 95% CI Lower Bound 95% CI Upper Bound    p Value
#>  2012.25 -0.02525882        -0.04586063       -0.004657013 0.04444444
#>  2012.50 -0.05117098        -0.07635097       -0.025990996 0.01111111
#>  2012.75 -0.04510912        -0.07028911       -0.019929129 0.01111111
#>  2013.00 -0.04426158        -0.06944157       -0.019081590 0.01111111
#>  2013.25 -0.05141523        -0.07659522       -0.026235237 0.01111111
#>  2013.50 -0.06877787        -0.09395786       -0.043597885 0.01111111
#>  2013.75 -0.05146447        -0.07664446       -0.026284478 0.01111111
#>  2014.00 -0.06947218        -0.09465217       -0.039714008 0.01111111
#>  2014.25 -0.06687073        -0.09662890       -0.037112564 0.01111111
#>  2014.50 -0.05345641        -0.08321457       -0.023698236 0.01111111
#>  2014.75 -0.04511139        -0.07486956       -0.015353223 0.02222222
#>  2015.00 -0.06350431        -0.09326248       -0.033746138 0.01111111
#>  2015.25 -0.05114785        -0.07632783       -0.025967857 0.01111111
#>  2015.50 -0.05904546        -0.08880363       -0.033865473 0.01111111
#>  2015.75 -0.05761259        -0.08737076       -0.027854422 0.01111111
#>  2016.00 -0.07352702        -0.10328519       -0.043768854 0.01111111
```


<img src="figure/fig_covsyn_resid-1.png" title="plot of chunk fig_covsyn_resid" alt="plot of chunk fig_covsyn_resid" style="display: block; margin: auto;" />


Finally, we can augment synth with many different outcome models. The simplest outcome model is a unit fixed effect model, which we can include by setting `fixedeff = T`.

```r

desyn <- augsynth(lngdpcapita ~ treated,
                   fips, year_qtr, kansas,
                   progfunc = "none", scm = T,
                   fixedeff = T)
```



```r
summary(desyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "none", scm = ..2, 
#>     fixedeff = ..3)
#> 
#> Average ATT Estimate (p Value for Joint Null):  -0.034   ( 0.328 )
#> L2 Imbalance: 0.082
#> Percent improvement from uniform weights: 55.1%
#> 
#> Avg Estimated Bias: -0.016
#> 
#> Inference type: Conformal inference
#> 
#>     Time    Estimate 95% CI Lower Bound 95% CI Upper Bound    p Value
#>  2012.25 -0.02151900        -0.04572669        0.005536651 0.07777778
#>  2012.50 -0.04567748        -0.06988516       -0.012925895 0.02222222
#>  2012.75 -0.03790080        -0.06210849       -0.005149216 0.04444444
#>  2013.00 -0.02357984        -0.04778753        0.003475817 0.07777778
#>  2013.25 -0.03260995        -0.05681764       -0.005554297 0.04444444
#>  2013.50 -0.05007896        -0.07428665       -0.023023308 0.02222222
#>  2013.75 -0.03457474        -0.05593447       -0.010367055 0.02222222
#>  2014.00 -0.04876841        -0.07297610       -0.018864793 0.02222222
#>  2014.25 -0.04661719        -0.07082488       -0.013865611 0.02222222
#>  2014.50 -0.03321017        -0.05741786        0.000000000 0.05555556
#>  2014.75 -0.02265951        -0.04686720        0.010092066 0.12222222
#>  2015.00 -0.03425329        -0.06130895        0.004194214 0.07777778
#>  2015.25 -0.02292051        -0.04712820        0.006983109 0.10000000
#>  2015.50 -0.02572370        -0.05277935        0.007027882 0.10000000
#>  2015.75 -0.02320911        -0.05026476        0.012390434 0.14444444
#>  2016.00 -0.03314605        -0.06589763        0.008149422 0.08888889
```


<img src="figure/fig_desyn-1.png" title="plot of chunk fig_desyn" alt="plot of chunk fig_desyn" style="display: block; margin: auto;" />

We can incorproate other outcome models by changing the `progfunc`.
Several outcome models are available, including, fitting the factor model directly with `gsynth`, general elastic net regression, bayesian structural time series estimation with `CausalImpact`, and matrix completion with `MCPanel`. For each outcome model you can supply an optional set of parameters, see documentation for details.


