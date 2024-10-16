## ----setup, include = FALSE------------------------------------------------------------------

knitr::opts_chunk$set(
    collapse = FALSE,
    message = FALSE,
    warning = FALSE,
    fig.align = 'center',
    comment = "#>"
)

## Install Synth if not already installed
# install.packages("Synth")

library(kableExtra)
library(magrittr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(rlang)

library(foreign)
library(graphics)
library(Synth)
library(augsynth)
data(kansas)

## For the Abadie plots
replication_theme <- theme_bw(base_family = "Times New Roman") + 
    theme(panel.grid = element_blank(), 
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0.0, 'pt'),
          legend.box.margin = margin(t = 5, r = 5, unit = 'pt'),
          legend.position = c(0.81, 0.9))

set.seed(1234)


## ----import tobacco data, include=FALSE------------------------------------------------------
# Get working data CA and 38 control states 
Wk.data <- read.dta(here::here("./tobacco_replication/data/smoking_wkdata_39.dta"))

# load in "Rest of US" data (computed from 50 states) for figure 1
RestUS.50 <- read.dta(here::here("./tobacco_replication/data/restofus_50.dta"))

# this data is used for table 1, where we average the 38 control states
RestUS.38 <- read.dta(here::here("./tobacco_replication/data/restofus_38.dta"))


## ----Synth dataprep and model, echo=FALSE, cache=TRUE, include=FALSE-------------------------
### Run Baseline Model

dataprep.out <- dataprep(
    foo = Wk.data,
    predictors = c("retprice", "xxincome", "K1_r_15_24"),
    predictors.op = c("mean"),
    dependent = c('cigsalepcap'),
    unit.variable = c('index'),
    time.variable = c('year'),
    special.predictors = list(list("xxbeer", 1984:1988, c("mean")),
                              list("cigsalepcap", 1988, c("mean")),
                              list("cigsalepcap", 1980, c("mean")),
                              list("cigsalepcap", 1975, c("mean"))
    ),
    treatment.identifier = 3,
    controls.identifier = c(1:39)[-3],
    time.predictors.prior = 1980:1988,
    time.optimize.ssr = 1970:1988,
    unit.names.variable = c('name'),
    time.plot = 1970:2005
)

### run synth

synth.out <- synth(
    dataprep.out,
    Pop.size = 300,
    Max.generations = 1,
    Unif.seed = 356987,
    Int.seed = 627478,
    Optim.method = "BFGS",
    L.ipop = 0,
    Maxiter.ipop = 1e10,
    Margin.ipop = 0.05,
    Sigf.ipop = 8,
    genoud = TRUE
)


## ----Abadie Figure 2, echo=FALSE, include=FALSE, fig.height=4, fig.width=4.6-----------------
f2_plot_df <- cbind(rownames(dataprep.out$Y1plot),
                    dataprep.out$Y1plot,
                    dataprep.out$Y0plot %*% synth.out$solution.w) %>%
    as_tibble() %>%
    mutate_all(as.numeric)

colnames(f2_plot_df) <- c('year', 'California', 'Synthetic California')

f2_plot_df <- f2_plot_df %>%
    reshape2::melt(id.vars = 'year',
                   variable.name = 'treat_group',
                   value.name = 'per_cap_cigs')

path_plot <- ggplot(f2_plot_df, aes(x = year, y = per_cap_cigs, linetype = treat_group)) +
    geom_line() +
    geom_vline(linetype = 'dotted', xintercept = 1988) +
    scale_linetype_manual(values = c('solid', 'dashed'),
                          breaks = c("California", "Synthetic California")) +
    scale_y_continuous(breaks = seq(0, 141, 20), limits = c(0, 140)) +
    scale_x_continuous(breaks = seq(1970, 2005, 5), limits = c(1970, 2005)) +
    labs(linetype = NULL,
         title = 'Abadie et al. (2010), Figure 2',
         y = 'per-capita cigarette sales (in packs)') +
    annotate('text', y = 41, x = 1987.5, label = 'Passage of Proposition 99 \u2192',
             hjust = 1, color = 'black', size = 3, family = "Times New Roman") +
    annotate('text', y = 0, x = 1970, label = 'ND replication: Synth package',
             hjust = 0, color = 'darkseagreen', size = 3) +
    replication_theme

path_plot


## ----Abadie Figure 3, echo=FALSE, include=FALSE, fig.height=4, fig.width=4.6-----------------
gap <- dataprep.out$Y1plot - dataprep.out$Y0plot %*% synth.out$solution.w
year <- dataprep.out$tag$time.plot

plot_df <- cbind(gap, year) %>% as_tibble()
colnames(plot_df) <- c('gap', 'year')

gap_plot <- ggplot(plot_df, aes(x = year, y = gap)) +
    geom_line() +
    geom_vline(linetype = 'dotted', xintercept = 1988) +
    geom_hline(linetype = 'dashed', yintercept = 0) +
    scale_y_continuous(breaks = seq(-30, 30, 10), limits = c(-30, 30)) +
    scale_x_continuous(breaks = seq(1970, 2005, 5), limits = c(1970, 2005)) +
    labs(linetype = NULL,
         title = 'Abadie et al. (2010), Figure 3',
         y = 'gap in per-capita cigarette sales (in packs)') +
    annotate('text', y = -25, x = 1987.5, label = 'Passage of Proposition 99 \u2192',
             hjust = 1, color = 'black', size = 3, family = "Times New Roman") +
    annotate('text', y = -30, x = 1970, label = 'ND replication: Synth package',
             hjust = 0, color = 'darkseagreen', size = 3) +
    replication_theme

gap_plot


## ----reformat tobacco data using tidyverse, echo=FALSE---------------------------------------
tobacco <- Wk.data %>%
    mutate(treated = ifelse((name == 'California') &
                                (year > 1988), 1, 0),
           state = name, # rename something meaningful
           state_index = index) %>%
    select(cigsalepcap, treated, state, year, retprice, xxincome, K1_r_15_24, xxbeer)

tobacco_70 <- tobacco %>% filter(year >= 1970) # subset data to 1970 forward to avoid issues due to missingess of outcomes


## ----run tobacco model using augsynth--------------------------------------------------------
# Basic SCM with covariates
syn <- augsynth(form = cigsalepcap ~ treated | retprice + xxincome + K1_r_15_24,
                unit = state,
                time = year,
                data = tobacco_70
                #progfunc = "none", scm=TRUE
)


## ----plot augsynth kansas models, echo=FALSE, cache=TRUE, fig.width=8, fig.height=4.5--------
all_results <- bind_rows(Conformal = summary(syn, inf_type = 'conformal')$att, 
                         Jackknife = summary(syn, inf_type = 'jackknife')$att, 
                         `Jackknife+` = summary(syn, inf_type = 'jackknife+')$att, 
                         Permutation = summary(syn, inf_type = 'permutation')$att, 
                         `Permutation (rstat)` = summary(syn, inf_type = 'permutation_rstat')$att, 
                         .id = 'name') %>% 
    mutate(name = factor(name, 
                         levels = c('Conformal', 'Jackknife', 'Jackknife+', "Permutation", "Permutation (rstat)"), 
                         ordered = T))

ggplot(all_results, aes(x = Time, y = Estimate, color = name)) + 
    geom_hline(yintercept = 0, linetype = 'solid') + 
    geom_vline(xintercept = syn$t_int, linetype = 'dashed') + 
    geom_line() +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = name), alpha = 0.2, size = 0.1) + 
    facet_wrap(. ~ name) + 
    scale_x_continuous(breaks = seq(1970, 2010, 5)) + 
    labs(x = 'Year', y = 'Gap in per-capita cigarette sales\n(in packs)', 
         caption = 'ND replication: augsynth package') +
    scale_color_manual(values = c('red', 'blue', 'darkgreen', 'purple', 'orange'), 
                       guide = 'none', aesthetics = c('color', 'fill')) +
    theme_bw() +
    theme(plot.caption = element_text(hjust = 0, color = 'darkseagreen', size = 9))


## ----show_perm_inference---------------------------------------------------------------------
syn_rstat = summary( syn, inf_type = "permutation_rstat" )
placebo_distribution( syn_rstat )


## ----fig.width=4.5, fig.height=3.25----------------------------------------------------------
plot(syn) + labs(title = 'Conformal inference (default)')


## ----fig.width=4.5, fig.height=4-------------------------------------------------------------
plot(syn_rstat, 
     plot_type = 'placebo',
     inf_type = 'permutation',
) + ggtitle("Replication of Abadie et al. (2010), Figure 4")


## ----fig.width=4.5, fig.height=3.25----------------------------------------------------------
p <- plot(syn_rstat, plot_type = 'estimate only') + ggtitle("Replication of Abadie et al. (2010), Figure 3")


## ----Format figures 3 comparison, echo=FALSE, fig.show="hold", fig.align='default', fig.height=4, fig.width=4.6----
gap_plot
p


## ----fig.width=4.5, fig.height=4-------------------------------------------------------------
p <- plot(syn_rstat, plot_type = 'outcomes') + 
    ggtitle("Replication of Abadie et al. (2010), Figure 2")


## ----Format figures 2 comparison, echo=FALSE, fig.show="hold", fig.align='default', fig.height=4, fig.width=4.6----
path_plot
p 


## ----Outcomes plot with average, echo=TRUE, fig.show="hold", fig.align='default', fig.height=4, fig.width=4.6----
plot(syn_rstat, plot_type="outcomes raw average")


## --------------------------------------------------------------------------------------------
donor_table(syn_rstat) %>%
    arrange( -abs(weight) )


## ----Figures 6 and 7, echo=TRUE, eval=TRUE, fig.show="hold", fig.align='default', fig.height=4.5, fig.width=4.6----
update_augsynth(syn) %>% # drops units with >20x treated RMSPE by default
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo' ) + 
    ggtitle("Replication of Abadie et al. (2010), Figure 5") + ylim(-51, 91) + 
    annotate('text', y = -48, x = 1970, label = "Removes donors with 20x California's  pre-treatment RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3) 
update_augsynth(syn, drop = 5) %>% 
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo') + 
    ggtitle("Replication of Abadie et al. (2010), Figure 6") + ylim(-51, 91) + 
    annotate('text', y = -48, x = 1970, label = "Removes donors with 5x California's pre-treatment RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3) 
update_augsynth(syn, drop = 2) %>% 
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo') + 
    ggtitle("Replication of Abadie et al. (2010), Figure 7") + ylim(-51, 91) + 
    annotate('text', y = -48, x = 1970, label = "Removes donors with 2x California's pre-treatment  RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3) 


## ----permutation without states, echo=TRUE, eval=TRUE, fig.show="hold", fig.align='default', fig.height=4.5, fig.width=4.6----
drop_states <- c("Iowa", "Arizona", "Alabama", "Illinois", "Indiana", "Idaho", "Connecticut", 
                 "New Mexico", "Texas", "Utah", "North Dakota", "South Dakota", "Vermont", 
                 "Wisconsin", "West Virginia", "Wyoming", "Tennessee", "Pennsylvania")

update_augsynth(syn, drop = drop_states) %>% 
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo') 

