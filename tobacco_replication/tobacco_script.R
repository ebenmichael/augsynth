
# Script replicating Abadie 2010 using augsynth

library( tidyverse )
library(augsynth)


set.seed(1234)


## ----reformat tobacco data using tidyverse, echo=FALSE-----------------------------
Wk.data <- read.dta(here::here("./tobacco_replication/data/smoking_wkdata_39.dta"))
tobacco <- Wk.data %>%
    mutate(treated = ifelse((name == 'California') &
                                (year > 1988), 1, 0),
           state = name, # rename something meaningful
           state_index = index) %>%
    select(cigsalepcap, treated, state, year, retprice, xxincome, K1_r_15_24, xxbeer)

tobacco_70 <- tobacco %>% filter(year >= 1970) # subset data to 1970 forward to avoid issues due to missingess of outcomes
rm( tobacco )
rm( Wk.data )

# Make beer variable
tobacco_70 <- tobacco_70 %>% group_by( state) %>%
    mutate( beer = mean( xxbeer[ year %in% 1984:1988 ], na.rm=TRUE ),
            retail = mean( retprice[ year %in% 1980:1988 ], na.rm=TRUE ),
            K1_r_15_24 = 100 * K1_r_15_24 ) %>%
    ungroup()


## ----run tobacco model using augsynth----------------------------------------------

rs = rnorm( length( unique( tobacco_70$state ) ), sd=0.2 )
names(rs) <- unique( tobacco_70$state )
tobacco_70$fnk = tobacco_70$year - 1970 + rs[ tobacco_70$state ]
summary( tobacco_70$fnk )
syn <- augsynth(form = cigsalepcap ~ treated | retail + xxincome + K1_r_15_24 + beer + fnk,
                unit = state,
                time = year,
                data = tobacco_70 )

ggplot( tobacco_70, aes( year, retprice, group = state ) ) +
    geom_line()
filter( tobacco_70, state == 'California' ) %>%
    dplyr::select( year, retprice ) %>%
    knitr::kable()
tobacco_70 %>%
    group_by( state ) %>%
    summarise( mn = mean( retprice ) )

covariate_balance_table( syn, pre_period = c( 1975, 1980, 1988 ) )

# Now inference
syn_rstat <- summary( syn, inf_type = 'permutation_rstat' )






## ----fig.width=4.5, fig.height=4---------------------------------------------------
p <- plot(syn_rstat, plot_type = 'outcomes raw average') +
    ggtitle("Replication of Abadie et al. (2010), Figure 2")
p


## ----------------------------------------------------------------------------------


## ----fig.width=4.5, fig.height=3.25------------------------------------------------
p <- plot(syn_rstat, plot_type = 'estimate only') +
    ggtitle("Replication of Abadie et al. (2010), Figure 3") +
    scale_y_continuous(limits = c(-30, 30))

p


## ----fig.width=4.5, fig.height=4---------------------------------------------------
plot(syn_rstat, plot_type = 'placebo' ) +
    ggtitle("Replication of Abadie et al. (2010), Figure 4")


## ----------------------------------------------------------------------------------
placebo_distribution( syn_rstat )


## ----------------------------------------------------------------------------------
donor_table(syn_rstat) %>%
    arrange( -abs(weight) )


## ----------------------------------------------------------------------------------
syn_classic <- augsynth(form = cigsalepcap ~ treated | retprice + xxincome + K1_r_15_24,
                unit = state,
                time = year,
                data = tobacco_70,
                progfunc = "none" )
donor_table(syn_classic) %>%
    filter( weight!= 0 ) %>%
    arrange( -abs(weight) )


## ----Figures 6 and 7, echo=TRUE, eval=TRUE, fig.show="hold", fig.align='default', fig.height=4.5, fig.width=4.6----
update_augsynth(syn) %>% # drops units with >20x treated RMSPE by default
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo' ) +
    ggtitle("Replication of Abadie et al. (2010), Figure 5") + ylim(-51, 91) +
    annotate('text', y = -48, x = 1970,
             label = "Removes donors with 20x California's  pre-treatment RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3)
update_augsynth(syn, drop = 5) %>%
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo') +
    ggtitle("Replication of Abadie et al. (2010), Figure 6") + ylim(-51, 91) +
    annotate('text', y = -48, x = 1970,
             label = "Removes donors with 5x California's pre-treatment RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3)
update_augsynth(syn, drop = 2) %>%
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo') +
    ggtitle("Replication of Abadie et al. (2010), Figure 7") + ylim(-51, 91) +
    annotate('text', y = -48, x = 1970,
             label = "Removes donors with 2x California's pre-treatment  RMSPE",
             hjust = 0, color = 'darkseagreen', size = 3)


## ----permutation without states, echo=TRUE, eval=TRUE, fig.show="hold", fig.align='default', fig.height=4.5, fig.width=4.6----
drop_states <- c("Iowa", "Arizona", "Alabama", "Illinois", "Indiana", "Idaho", "Connecticut",
                 "New Mexico", "Texas", "Utah", "North Dakota", "South Dakota", "Vermont",
                 "Wisconsin", "West Virginia", "Wyoming", "Tennessee", "Pennsylvania")

update_augsynth(syn, drop = drop_states) %>%
    summary( inf_type = 'permutation' ) %>%
    plot(plot_type = 'placebo')


rm( p )
