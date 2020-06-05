#' Economic indicators for US states from 1990-2016
#' 
#' 
#' @format A dataframe with 5250 rows and 32 variables:
#' \describe{
#'  \item{fips}{FIPS code for each state}
#'  \item{year}{Year of measurement}
#'  \item{qtr}{Quarter (1-4) of measurement}
#'  \item{state}{Name of State}
#'  \item{gdp}{Gross State Product (millions of $) Values before 2005 are linearly interpolated between years}
#'  \item{revenuepop}{State and local revenue per capita}
#'  \item{rev_state_total}{State total general revenue (millions of $)}
#'  \item{rev_local_total}{Local total general revenue (millions of $)}
#'  \item{popestimate}{Population estimate}
#'  \item{qtrly_estabs_count}{Count of establishments for a given quarter}
#'  \item{month1_emplvl, month2_emplvl, month3_emplvl}{ Employment level for first, second, and third months of a given quarter}
#'  \item{total_qtrly_wages}{Total wages for a givne quarter}
#'  \item{taxable_qtrly_wage}{Taxable wages for a given quarter}
#'  \item{avg_wkly_wage}{Average weekly wage for a given quarter}
#'  \item{year_qtr}{Year and quarter combined into one continuous variable}
#'  \item{treated}{Whether the state passed tax cuts before the given year and quareter}
#'  \item{lngdpcapita}{Natural log of GDP per capita}
#'  \item{emplvlcapita}{Average employment level per capita}
#'  \item{Xcapita}{Per capita value of X}
#'  \item{abb}{State abbreviation}
#' }
"kansas"