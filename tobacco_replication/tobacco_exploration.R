

# Exploration of the tobacco data


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


## For the Abadie plots
replication_theme <- theme_bw(base_family = "Times New Roman") +
    theme(panel.grid = element_blank(),
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0.0, 'pt'),
          legend.box.margin = margin(t = 5, r = 5, unit = 'pt'),
          legend.position = c(0.81, 0.9))

set.seed(1234)


# Get working data CA and 38 control states
Wk.data <- read.dta(here::here("./tobacco_replication/data/smoking_wkdata_39.dta"))



time_tx = 1988


# Look at the data ----

skimr::skim( Wk.data )

# How much do the predictors predict our outcome?
library( lme4 )
Wk.data <- mutate( Wk.data,
                   year_c = year - 1988 )

M <- lmer( cigsalepcap ~ 1 + year_c + I(year_c^2) + retprice + xxincome + K1_r_15_24 + (1|index),
           data = Wk.data )

arm::display( M )


d1975 = filter( Wk.data, year==1975 )
table( d1975$index )

filter( Wk.data, year == 1975 ) %>%
    dplyr::select( cigsalepcap ) %>%
    head()

# name and index correspond - check
table( table( Wk.data$name, Wk.data$index ) )

al = filter( Wk.data, name == "Alabama")
al %>%
    summarise( xxbeer = mean( xxbeer, na.rm=TRUE ),
               rp = sd( retprice, na.rm=TRUE ),
               xx = sd( xxincome, na.rm=TRUE ),
               #cig1989 = cigsalepcap[ year == 1989 ],
               cig1980 = length( cigsalepcap[ year == 1980 ] ),
               cig1988 = length( cigsalepcap[ year == 1988 ] ),
               cig1975 = length( cigsalepcap[ year == 1975 ] ),
               .groups = "drop" )


covs <- Wk.data %>%
    #dplyr::filter( year %in% 1984:1988 ) %>%
    group_by( name, index ) %>%
    summarise( xxbeer = mean( xxbeer[ year %in% 1984:1988 ], na.rm=TRUE ),
               rp = sd( retprice, na.rm=TRUE ),
               xx = sd( xxincome, na.rm=TRUE ),
               retprice = mean( retprice, na.rm=TRUE ),
               xxincome = mean( xxincome, na.rm=TRUE ),
               K1_r_15_24 = mean( K1_r_15_24, na.rm=TRUE ),
               cig1989 = cigsalepcap[ year == 1989 ],
               cig1995 = cigsalepcap[ year == 1989 ],
               cig1980 =  cigsalepcap[ year == 1980 ] ,

               cig1988 =  cigsalepcap[ year == 1988 ] ,
               cig1975 =  cigsalepcap[ year == 1975 ] ,
               .groups = "drop"
    )
head( covs )



MRpre <- lm( cig1995 ~ cig1980 + cig1988 + cig1975, data=covs )
MRcov <- lm( cig1995 ~ xxbeer + retprice + xxincome + K1_r_15_24, data=covs )

summary( MRpre )
summary( MRcov )

MRfull <- lm( cig1995 ~ cig1980 + cig1988 + cig1975 + xxbeer + retprice + xxincome + K1_r_15_24, data=covs )
summary( MRfull )

anova( MRpre, MRfull )



## ----reformat tobacco data using tidyverse ---------------------------------------
tobacco <- Wk.data %>%
    mutate(treated = ifelse((name == 'California') &
                                (year > 1988), 1, 0),
           state = name, # rename something meaningful
           state_index = index) %>%
    select(cigsalepcap, treated, state, year, retprice, xxincome, K1_r_15_24, xxbeer)


# subset data to 1970 forward to avoid issues due to missingess of outcomes
tobacco_70 <- tobacco %>% filter(year >= 1970)


tobacco %>% group_by( year ) %>%
    summarise( prop.miss = mean( is.na( cigsalepcap ))) %>%
    print( n=100 )

tobacco_55 <- tobacco %>% filter(year >= 1955)


# Plot of states outcome over time ----

ggplot( Wk.data, aes( year, cigsalepcap, group=index ) ) +
    geom_line() +
    scale_y_log10()


# Various plots of covariates over time ----
ggplot( Wk.data, aes( year, K1_r_15_24, group=index ) ) +
    geom_line()

ggplot( Wk.data, aes( year, xxincome, group=index ) ) +
    geom_line()

ggplot( Wk.data, aes( year, retprice, group=index ) ) +
    geom_line() +
    scale_y_log10()

ggplot( Wk.data, aes( year, cigsalepcap, group=index ) ) +
    geom_line() +
    scale_y_log10()

ggplot( Wk.data, aes( year, xxbeer, group=index ) ) +
    geom_line() +
    facet_wrap( ~ name ) +
    scale_y_log10()

ggplot( Wk.data, aes( xxbeer ) ) +
    geom_histogram()


# Windsorize the units to remove extreme outliers


windsorize_outliers <- function(column) {
    Q1 <- quantile(column, 0.25, na.rm = TRUE)
    Q3 <- quantile(column, 0.75, na.rm = TRUE)
    IQR_value <- IQR(column, na.rm = TRUE)

    lower_bound <- Q1 - 2 * IQR_value
    upper_bound <- Q3 + 2 * IQR_value

    # Windsorize outliers
    column <- ifelse(column < lower_bound, lower_bound,
                     ifelse(column > upper_bound, upper_bound, column))
    return(column)
}

replace_outliers_with_na <- function(column) {
    Q1 <- quantile(column, 0.25, na.rm = TRUE)
    Q3 <- quantile(column, 0.75, na.rm = TRUE)
    IQR_value <- IQR(column, na.rm = TRUE)

    lower_bound <- Q1 - 2 * IQR_value
    upper_bound <- Q3 + 2 * IQR_value

    prop = column < lower_bound | column > upper_bound
    cat( glue::glue( "{deparse(substitute(column))} outliers %" ),
         100*mean(prop, na.rm=TRUE), "\n" )
    column[column < lower_bound | column > upper_bound] <- NA
    return(column)
}


# Apply the outlier function to each specified column
data = Wk.data
data <- data %>%
    mutate(across(c(xxincome, xxbeer, K1_r_15_24,
                    rate, cigsale, cigsalepcap, retprice,
                    pop_census), replace_outliers_with_na))


ggplot( Wk.data, aes( xxbeer ) ) +
    geom_histogram()

ggplot( data, aes( xxbeer ) ) +
    geom_histogram()


ggplot( Wk.data, aes( rate ) ) +
    geom_histogram()

ggplot( data, aes( rate ) ) +
    geom_histogram()





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

if ( FALSE ) {
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
    saveRDS( synth.out, file = here::here( "tobacco_replication/Synth_result.rds" ) )

}

synth.out = readRDS( here::here( "tobacco_replication/Synth_result.rds" ) )


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



## ----run tobacco model using augsynth--------------------------------------------------------

# Basic SCM with covariates
syn <- augsynth(form = cigsalepcap ~ treated | retprice + xxincome + K1_r_15_24,
                unit = state,
                time = year,
                data = tobacco_70 )

sum = summary(syn, inf_type = 'permutation_rstat')


# Look at analysis with shifting initial year of treatment ----

do_year <- function( tob, yr ) {
    tob <- tob %>%
        dplyr::mutate( tx_yr = ifelse((state == 'California') & (year > yr), 1, 0) )

    sum2 <- augsynth(form = cigsalepcap ~ tx_yr | retprice + xxincome + K1_r_15_24,
                     unit = state,
                     time = year,
                     data = tob ) %>%
        summary(inf_type = 'permutation_rstat')

    sum2$att$tx_year = yr
    sum2$att
}

library( tidyverse )
yrs <- map_df( 1975:1988, do_year, tob = tobacco_55 )
#yrs$tx_year = as.factor(yrs$tx_year)
head( yrs )

ggplot( yrs, aes(x = Time, y = Estimate, color = tx_year, group=tx_year)) +
    geom_hline(yintercept = 0, linetype = 'solid') +
    geom_vline(xintercept = syn$t_int - 0.5, linetype = 'dashed') +
    geom_line() +
#    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = name), alpha = 0.2, size = 0.1) +
    #facet_wrap( ~ name) +
    scale_x_continuous(breaks = seq(1970, 2010, 5)) +
    labs(x = 'Year', y = 'Gap in per-capita cigarette sales\n(in packs)' ) +
    theme_bw() +
    theme(plot.caption = element_text(hjust = 0, color = 'darkseagreen', size = 9))







ggplot(all_results, aes(x = Time, y = Estimate, color = name)) +
    geom_hline(yintercept = 0, linetype = 'solid') +
    geom_vline(xintercept = syn$t_int - 0.5, linetype = 'dashed') +
    geom_line() +
    #geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound, fill = name), alpha = 0.2, size = 0.1) +
    #facet_wrap( ~ name) +
    scale_x_continuous(breaks = seq(1970, 2010, 5)) +
    labs(x = 'Year', y = 'Gap in per-capita cigarette sales\n(in packs)' ) +
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

