library(haven)
library(tidyverse)

kansas <- read_dta("kansas_longer2.dta")
state_abb <- read_csv("us-state-ansi-fips.csv") %>%
                rename(fips = st, abb = stusps) %>%
                mutate(fips = as.numeric(fips)) %>%
                select(fips, abb)
kansas <- kansas %>%
        rename(fips=Fips) %>%
        filter(year >= 1990,
               !is.na(fips), # filter out all of US
               fips != 11, # filter out DC
            #    year_qtr >= 2005 | year_qtr == round(year_qtr)
               ) %>%
        # interpolate GDP
        mutate(year_qtr = year + qtr / 4 - 0.25, # combine year and quarter
               fips = as.integer(fips), # state id
               treated = 1 * (fips == 20) * (year_qtr >= 2012.25),
               gdp = ifelse((qtr == 1) | (year >= 2005), gdp, NA),
               popestimate = ifelse((qtr == 1), popestimate, NA)) %>%
        # interpolate GDP and population
        group_by(fips) %>%
        arrange(year_qtr) %>%
        mutate(gdp = approx(year_qtr, gdp, year_qtr)$y,
               popestimate = approx(year_qtr, popestimate, year_qtr)$y) %>%
        ungroup() %>% arrange(fips, year_qtr) %>%
        mutate(gdpcapita = gdp / popestimate * 1e6,
               lngdp = log(gdp),
               lngdpcapita = log(gdpcapita),
               revstatecapita = rev_state_total / popestimate * 1e6,
               revlocalcapita = rev_local_total / popestimate * 1e6,
               emplvl1capita = month1_emplvl / popestimate,
               emplvl2capita = month2_emplvl / popestimate,
               emplvl3capita = month3_emplvl / popestimate,
               emplvlcapita = (month1_emplvl + month2_emplvl + month3_emplvl) / (3 * popestimate),
               totalwagescapita = total_qtrly_wages / popestimate,
               taxwagescapita = taxable_qtrly_wages / popestimate,
               avgwklywagecapita = avg_wkly_wage,
               estabscapita = qtrly_estabs_count / popestimate) %>%
        filter(year_qtr <= 2016) %>%
        inner_join(state_abb)

for (name in colnames(kansas)) {
        attributes(kansas[[name]])$label = NULL 
}

