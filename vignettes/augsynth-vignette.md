---
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
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


## Example: The economic costs of conflict

### The data
To show the usage and features of `augsynth`, we'll use data on the impact of terrorism in the Basque country supplied in the `Synth` package. Our interest is in estimating the effect of increased terrorism in the Basque country on GDP per capita.



```r
library(magrittr)
library(dplyr)
library(Synth)
data(basque)
```

The `basque` dataset contains the GDP per capita (the outcome measure) for the regions of Spain from 1955 to 1997, as well as data on several auxiliary covariates at a few time points in the study period.
<table class="table table-hover table-responsive" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> regionno </th>
   <th style="text-align:left;"> regionname </th>
   <th style="text-align:right;"> year </th>
   <th style="text-align:right;"> gdpcap </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1955 </td>
   <td style="text-align:right;"> 2.354542 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1956 </td>
   <td style="text-align:right;"> 2.480149 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1957 </td>
   <td style="text-align:right;"> 2.603613 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1958 </td>
   <td style="text-align:right;"> 2.637104 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1959 </td>
   <td style="text-align:right;"> 2.669880 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Spain (Espana) </td>
   <td style="text-align:right;"> 1960 </td>
   <td style="text-align:right;"> 2.869966 </td>
  </tr>
</tbody>
</table>


To run `augsynth`, we need to include a treatment status column that indicates which region was treated and at what time. As in the original study, we'll mark 1975 as the start of increased terrorism in the Basque country. We'll also drop the data for the whole country of Spain.


```r
basque %>%
    mutate(
        treatment=case_when(year < 1975 ~ 0,
                            regionno != 17 ~ 0,
                            regionno == 17 ~ 1) # Basque after 1975 is treated
    ) %>%
    filter(regionno != 1) -> basque
```


### Synth
Now to find a synthetic control using the entire series of pre-intervention outcomes (and no auxiliary covariates), we can use `augsynth`. To do so we just need to give `augsynth` a formula like `outcome ~ treatment`, tell it what the unit and time variables are, optionally provide when intervention took place (the code will automatically determine this if `t_int` is not provided), and specify that we don't want to fit an outcome model


```r
library(augsynth)
syn <- augsynth(gdpcap ~ treatment, regionno, year, basque,
                progfunc="None", scm=T, t_int=1975)
```

We can then look at the ATT estimates for each post-intervention time period and overall. We'll also see standard errors estimated using leave-out-one estimates of the noise and the quality of the synthetic control fit measured by the L2 distance between Basque and its synthetic control.


```r
summary(syn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "None", scm = ..2)
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.674  (0.108)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.377  (0.052)	Avg Estimated Bias: NA
#> 
#>  Time    Estimate  Std.Error
#>  1975  0.14624071 0.09479813
#>  1976  0.01401143 0.07851549
#>  1977 -0.11450644 0.06315419
#>  1978 -0.28092436 0.08247284
#>  1979 -0.40972206 0.06435448
#>  1980 -0.56006344 0.08661734
#>  1981 -0.72817478 0.09819208
#>  1982 -0.76861559 0.08676485
#>  1983 -0.80863120 0.08030836
#>  1984 -0.73939424 0.14511637
#>  1985 -0.66004014 0.21879683
#>  1986 -0.76871587 0.20744798
#>  1987 -0.85959322 0.19475415
#>  1988 -0.93892210 0.13390084
#>  1989 -1.00259722 0.08098128
#>  1990 -0.97145978 0.05072853
#>  1991 -0.96544129 0.02748016
#>  1992 -0.92196338 0.06679979
#>  1993 -0.88307019 0.11571515
#>  1994 -0.92080987 0.04759818
#>  1995 -0.80734361 0.03887265
#>  1996 -0.79383498 0.04024409
#>  1997 -0.76478286 0.07001232
```

It's easier to see this information visually. Below we plot the difference between the Basque region and it's synthetic control. Before the increase in terrorism (to the left of the dashed line) we expect these to be close, and after the increase we measure the effect (plus or minus 2 standard errors).


```r
plot(syn)
```

<img src="figure/fig_syn-1.png" title="plot of chunk fig_syn" alt="plot of chunk fig_syn" style="display: block; margin: auto;" />

### Augmenting synth with an outcome model
In this example the pre-intervention synthetic control fit is quite good: the L2 imbalance is 0.377, about 5% of the imbalance between the Basque country and the average of the other regions. We can get slightly by _augmenting_ synth with ridge regression. To do this we change `progfunc` to `"Ridge"`. We can also choose the ridge hyper-parameter by setting `lambda`:

```r
asyn <- augsynth(gdpcap ~ treatment, regionno, year, basque,
                progfunc="Ridge", scm=T, lambda=8)
```

We can look at the summary and plot the results. Now in the summary output we see an estimate of the overall bias of synth; we measure this with the average amount that augmentation changes the synth estimate. Notice that the estimates don't change very much, but the standard errors are tighter.

```r
summary(asyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "Ridge", scm = ..2, 
#>     lambda = 8)
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.680  (0.177)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.373  (0.051)	Avg Estimated Bias: 0.006
#> 
#>  Time    Estimate  Std.Error
#>  1975  0.14284671 0.03043658
#>  1976  0.01051501 0.04488431
#>  1977 -0.11818961 0.05996912
#>  1978 -0.28636354 0.08893771
#>  1979 -0.41587619 0.08393560
#>  1980 -0.56561235 0.09575623
#>  1981 -0.73282714 0.10428851
#>  1982 -0.77309648 0.11628085
#>  1983 -0.81308151 0.13352415
#>  1984 -0.74408923 0.18093616
#>  1985 -0.66512822 0.25148433
#>  1986 -0.77262185 0.24281910
#>  1987 -0.86267153 0.23462411
#>  1988 -0.94170897 0.19053099
#>  1989 -1.00523288 0.16744572
#>  1990 -0.97578145 0.17977609
#>  1991 -0.97126844 0.19776985
#>  1992 -0.92882622 0.19679309
#>  1993 -0.89091768 0.21206961
#>  1994 -0.92981430 0.21414406
#>  1995 -0.81770328 0.23041534
#>  1996 -0.80406052 0.23205106
#>  1997 -0.77580999 0.25197426
```


```r
plot(asyn)
```


<img src="figure/fig_asyn-1.png" title="plot of chunk fig_asyn" alt="plot of chunk fig_asyn" style="display: block; margin: auto;" />



There are also several auxiliary covariates. We can include these in the augmentation by fitting an outcome model using the auxiliary covariates. To do this we simply add the covariates into the formula after `|`; by default this will average the auxiliary covariates over the pre-intervention period, dropping `NA` values and regress out the auxiliary covariates.

```r
covsyn <- augsynth(gdpcap ~ treatment | invest + sec.agriculture + sec.energy + gdpcap,
                   regionno, year, basque,
                   progfunc="None", scm=T)
```

Again we can look at the summary and plot the results.

```r
summary(covsyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "None", scm = ..2)
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.598  (0.584)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.356  (0.049)	Avg Estimated Bias: NA
#> 
#>  Time    Estimate Std.Error
#>  1975  0.08079169 0.1694155
#>  1976 -0.09032869 0.1741993
#>  1977 -0.25127091 0.1852070
#>  1978 -0.50497558 0.1287235
#>  1979 -0.63415654 0.1065111
#>  1980 -0.75488930 0.1345742
#>  1981 -0.88888610 0.2304046
#>  1982 -0.87390861 0.3200002
#>  1983 -0.87595638 0.3968439
#>  1984 -0.83598023 0.3929884
#>  1985 -0.80745513 0.3769219
#>  1986 -0.79428790 0.5110267
#>  1987 -0.77457014 0.6328869
#>  1988 -0.75849526 0.7372358
#>  1989 -0.73422543 0.8302636
#>  1990 -0.64795748 0.8650987
#>  1991 -0.60181330 0.8930641
#>  1992 -0.60586198 0.7674044
#>  1993 -0.61433586 0.6476947
#>  1994 -0.62430275 0.7368969
#>  1995 -0.49413970 0.8024834
#>  1996 -0.37775035 0.8398268
#>  1997 -0.29207086 0.8760892
```


```r
plot(covsyn)
```

<img src="figure/fig_covsyn-1.png" title="plot of chunk fig_covsyn" alt="plot of chunk fig_covsyn" style="display: block; margin: auto;" />

Now we can additionally fit ridge ASCM on the residuals, look at the summary, and plot the results.

```r
covsyn_aug <- augsynth(gdpcap ~ treatment | invest + sec.agriculture + sec.energy + gdpcap,
                   regionno, year, basque,
                   progfunc="Ridge", scm=T, lambda = 1e-1)
```


```r
summary(covsyn_aug)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "Ridge", scm = ..2, 
#>     lambda = 0.1)
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.539  (0.594)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.288  (0.040)	Avg Estimated Bias: 1.482
#> 
#>  Time    Estimate  Std.Error
#>  1975  0.08995881 0.12205497
#>  1976 -0.08746053 0.12894165
#>  1977 -0.25656767 0.14431761
#>  1978 -0.50864586 0.12886892
#>  1979 -0.64099888 0.07574878
#>  1980 -0.74419804 0.06716698
#>  1981 -0.84795886 0.18677276
#>  1982 -0.81472171 0.27522913
#>  1983 -0.80217348 0.35103027
#>  1984 -0.75020946 0.30219350
#>  1985 -0.71104319 0.24401438
#>  1986 -0.66245109 0.38183020
#>  1987 -0.61625347 0.51655679
#>  1988 -0.60298344 0.65509943
#>  1989 -0.58387592 0.78779500
#>  1990 -0.53007566 0.86164336
#>  1991 -0.51419468 0.91339218
#>  1992 -0.54541279 0.79708592
#>  1993 -0.58157345 0.68305437
#>  1994 -0.57517346 0.81015433
#>  1995 -0.46677923 0.90654202
#>  1996 -0.35795188 0.97378675
#>  1997 -0.28411469 1.03518310
```


```r
plot(covsyn_aug)
```

<img src="figure/fig_covsyn_aug-1.png" title="plot of chunk fig_covsyn_aug" alt="plot of chunk fig_covsyn_aug" style="display: block; margin: auto;" />


Finally, we can augment synth with many different outcome models, this is as easy as changing the `progfunc`. For instance, we can augment synth with matrix completion using `mcpanel`.


```r
mcpsyn <- augsynth(gdpcap ~ treatment,
                   regionno, year, basque,
                   progfunc="MCP", scm=T)
```

For the other outcome models we do not (yet) supply standard error estimates.

```r
summary(mcpsyn)
#> 
#> Call:
#> single_augsynth(form = form, unit = !!enquo(unit), time = !!enquo(time), 
#>     t_int = t_int, data = data, progfunc = "MCP", scm = ..2)
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.614  (NA)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.377  (0.052)	Avg Estimated Bias: -0.032
#> 
#>  Time    Estimate Std.Error
#>  1975  0.10108009        NA
#>  1976 -0.05014203        NA
#>  1977 -0.20185743        NA
#>  1978 -0.48659799        NA
#>  1979 -0.64054449        NA
#>  1980 -0.75309182        NA
#>  1981 -0.83800892        NA
#>  1982 -0.85905119        NA
#>  1983 -0.87544758        NA
#>  1984 -0.85608496        NA
#>  1985 -0.83126554        NA
#>  1986 -0.81854095        NA
#>  1987 -0.79357362        NA
#>  1988 -0.78314084        NA
#>  1989 -0.75117983        NA
#>  1990 -0.67507896        NA
#>  1991 -0.62609311        NA
#>  1992 -0.67822034        NA
#>  1993 -0.72995631        NA
#>  1994 -0.67069581        NA
#>  1995 -0.51743996        NA
#>  1996 -0.43540996        NA
#>  1997 -0.36030222        NA
```


```r
plot(mcpsyn)
```

<img src="figure/fig_mcpsyn-1.png" title="plot of chunk fig_mcpsyn" alt="plot of chunk fig_mcpsyn" style="display: block; margin: auto;" />

Several other outcome models are available, including general elastic net regression, bayesian structural time series estimation with `CausalImpact`, and the generalized synthetic control method `gsynth`. For each outcome model you can supply an optional set of parameters, see documentation for details.
