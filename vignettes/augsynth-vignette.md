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
Now to find a synthetic control using the entire series of pre-intervention outcomes (and no auxiliary covariates), we can use `augsynth`. To do so we just need to give `augsynth` a formula like `outcome ~ treatment`, tell it what the unit and time variables are and when the intervention took place, and specify that we don't want to fit an outcome model


```r
library(augsynth)
syn <- augsynth(gdpcap ~ treatment, regionno, year, 1975, basque,
                progfunc="None", weightfunc="SCM")
```

We can then look at the ATT estimates for each post-intervention time period and overall. We'll also see standard errors estimated using leave-out-one estimates of the noise and the quality of the synthetic control fit measured by the L2 distance between Basque and its synthetic control.


```r
summary(syn)
#> 
#> Call:
#> augsynth(form = gdpcap ~ treatment, unit = regionno, time = year, 
#>     t_int = 1975, data = basque, progfunc = "None", weightfunc = "SCM")
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.674  (0.769)
#> Std. Deviation: 0.594
#> 
#> L2 Imbalance (Scaled): 0.377  (0.052)	Avg Estimated Bias: NA
#> 
#>  Time    Estimate Std.Error
#>  1975  0.14624071 0.5009803
#>  1976  0.01401143 0.5161356
#>  1977 -0.11450644 0.5439447
#>  1978 -0.28092436 0.5246125
#>  1979 -0.40972206 0.4850378
#>  1980 -0.56006344 0.5082463
#>  1981 -0.72817478 0.5441567
#>  1982 -0.76861559 0.5705742
#>  1983 -0.80863120 0.6136181
#>  1984 -0.73939424 0.6946066
#>  1985 -0.66004014 0.7904961
#>  1986 -0.76871587 0.8364455
#>  1987 -0.85959322 0.8805142
#>  1988 -0.93892210 0.8755663
#>  1989 -1.00259722 0.8855021
#>  1990 -0.97145978 0.9105807
#>  1991 -0.96544129 0.9594822
#>  1992 -0.92196338 0.9060451
#>  1993 -0.88307019 0.8677904
#>  1994 -0.92080987 0.8929513
#>  1995 -0.80734361 0.9250645
#>  1996 -0.79383498 0.9620763
#>  1997 -0.76478286 0.9955493
```

It's easier to see this information visually. Below we plot the difference between the Basque region and it's synthetic control. Before the increase in terrorism (to the left of the dashed line) we expect these to be close, and after the increase we measure the effect (plus or minus 2 standard errors).


```r
plot(syn)
```

<img src="figure/fig_syn-1.png" title="plot of chunk fig_syn" alt="plot of chunk fig_syn" style="display: block; margin: auto;" />

### Augmenting synth with an outcome model
In this example the pre-intervention synthetic control fit is quite good: the L2 imbalance is 0.377, about 5% of the imbalance between the Basque country and the average of the other regions. We can get slightly by _augmenting_ synth with ridge regression. To do this we change `progfunc` to `"Ridge"`. We can also choose the ridge hyper-parameter by setting `lambda` in the option list of outcome model arguments `opts_out`:

```r
asyn <- augsynth(gdpcap ~ treatment, regionno, year, 1975, basque,
                progfunc="Ridge", weightfunc="SCM", opts_out=list(lambda=8))
```

We can look at the summary and plot the results. Now in the summary output we see an estimate of the overall bias of synth; we measure this with the average amount that augmentation changes the synth estimate. Notice that the estimates don't change very much, but the standard errors are tighter.

```r
summary(asyn)
#> 
#> Call:
#> augsynth(form = gdpcap ~ treatment, unit = regionno, time = year, 
#>     t_int = 1975, data = basque, progfunc = "Ridge", weightfunc = "SCM", 
#>     opts_out = list(lambda = 8))
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.680  (0.64)
#> Std. Deviation: 0.495
#> 
#> L2 Imbalance (Scaled): 0.373  (0.051)	Avg Estimated Bias: 0.006
#> 
#>  Time    Estimate Std.Error
#>  1975  0.14284671 0.2547758
#>  1976  0.01051501 0.2989904
#>  1977 -0.11818961 0.3555713
#>  1978 -0.28636354 0.3440486
#>  1979 -0.41587619 0.3066112
#>  1980 -0.56561235 0.3223258
#>  1981 -0.73282714 0.3687233
#>  1982 -0.77309648 0.3983652
#>  1983 -0.81308151 0.4487968
#>  1984 -0.74408923 0.5045262
#>  1985 -0.66512822 0.5821369
#>  1986 -0.77262185 0.6380248
#>  1987 -0.86267153 0.6934532
#>  1988 -0.94170897 0.6897501
#>  1989 -1.00523288 0.7098273
#>  1990 -0.97578145 0.7679853
#>  1991 -0.97126844 0.8464035
#>  1992 -0.92882622 0.8014459
#>  1993 -0.89091768 0.7718662
#>  1994 -0.92981430 0.8285423
#>  1995 -0.81770328 0.8760488
#>  1996 -0.80406052 0.9364995
#>  1997 -0.77580999 0.9968989
```


```r
plot(asyn)
```


<img src="figure/fig_asyn-1.png" title="plot of chunk fig_asyn" alt="plot of chunk fig_asyn" style="display: block; margin: auto;" />



There are also several auxiliary covariates. We can include these in the augmentation by fitting an outcome model using the auxiliary covariates. To do this we simply add the covariates into the formula after `|`; by default this will average the auxiliary covariates over the pre-intervention period, dropping `NA` values. We also set `lambda=0` to fit OLS.

```r
covsyn <- augsynth(gdpcap ~ treatment | invest + sec.agriculture + sec.energy + gdpcap,
                   regionno, year, 1975, basque,
                   progfunc="Ridge", weightfunc="SCM", opts_out=list(lambda=0))
```

Again we can look at the summary and plot the results.

```r
summary(covsyn)
#> 
#> Call:
#> augsynth(form = gdpcap ~ treatment | invest + sec.agriculture + 
#>     sec.energy + gdpcap, unit = regionno, time = year, t_int = 1975, 
#>     data = basque, progfunc = "Ridge", weightfunc = "SCM", opts_out = list(lambda = 0))
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.598  (0.673)
#> Std. Deviation: 0.522
#> 
#> L2 Imbalance (Scaled): 0.356  (0.049)	Avg Estimated Bias: 1.541
#> 
#>  Time    Estimate Std.Error
#>  1975  0.08079169 0.2794998
#>  1976 -0.09032869 0.3171225
#>  1977 -0.25127091 0.3644139
#>  1978 -0.50497558 0.3327159
#>  1979 -0.63415654 0.2724416
#>  1980 -0.75488930 0.2659063
#>  1981 -0.88888610 0.3278893
#>  1982 -0.87390861 0.3980760
#>  1983 -0.87595638 0.4833156
#>  1984 -0.83598023 0.5183797
#>  1985 -0.80745513 0.5608427
#>  1986 -0.79428790 0.6477533
#>  1987 -0.77457014 0.7305508
#>  1988 -0.75849526 0.7806824
#>  1989 -0.73422543 0.8401228
#>  1990 -0.64795748 0.8573010
#>  1991 -0.60181330 0.8961347
#>  1992 -0.60586198 0.8462190
#>  1993 -0.61433586 0.8189529
#>  1994 -0.62430275 0.9037494
#>  1995 -0.49413970 0.9626462
#>  1996 -0.37775035 0.9524442
#>  1997 -0.29207086 0.9638151
```


```r
plot(covsyn)
```

<img src="figure/fig_covsyn-1.png" title="plot of chunk fig_covsyn" alt="plot of chunk fig_covsyn" style="display: block; margin: auto;" />


Finally, we can augment synth with many different outcome models, this is as easy as changing the `progfunc`. For instance, we can augment synth with matrix completion using `mcpanel`.


```r
mcpsyn <- augsynth(gdpcap ~ treatment,
                   regionno, year, 1975, basque,
                   progfunc="MCP", weightfunc="SCM")
```

For the other outcome models we do not (yet) supply standard error estimates.

```r
summary(mcpsyn)
#> 
#> Call:
#> augsynth(form = gdpcap ~ treatment, unit = regionno, time = year, 
#>     t_int = 1975, data = basque, progfunc = "MCP", weightfunc = "SCM")
#> 
#> Average ATT Estimate (Pooled Std. Error): -0.614  (NA)
#> Std. Deviation: NA
#> 
#> L2 Imbalance (Scaled): 0.377  (0.052)	Avg Estimated Bias: -0.060
#> 
#>  Time   Estimate Std.Error
#>  1975  0.1010658        NA
#>  1976 -0.0501578        NA
#>  1977 -0.2018750        NA
#>  1978 -0.4866078        NA
#>  1979 -0.6405479        NA
#>  1980 -0.7530885        NA
#>  1981 -0.8379952        NA
#>  1982 -0.8590414        NA
#>  1983 -0.8754413        NA
#>  1984 -0.8560811        NA
#>  1985 -0.8312636        NA
#>  1986 -0.8185365        NA
#>  1987 -0.7935674        NA
#>  1988 -0.7831366        NA
#>  1989 -0.7511763        NA
#>  1990 -0.6750793        NA
#>  1991 -0.6260993        NA
#>  1992 -0.6782318        NA
#>  1993 -0.7299734        NA
#>  1994 -0.6707021        NA
#>  1995 -0.5174496        NA
#>  1996 -0.4354104        NA
#>  1997 -0.3603028        NA
```


```r
plot(mcpsyn)
```

<img src="figure/fig_mcpsyn-1.png" title="plot of chunk fig_mcpsyn" alt="plot of chunk fig_mcpsyn" style="display: block; margin: auto;" />

Several other outcome models are available, including general elastic net regression, bayesian structural time series estimation with `CausalImpact`, and the generalized synthetic control method `gsynth`. For each outcome model you can supply an optional set of parameters `opts_out`, see documentation for details.
