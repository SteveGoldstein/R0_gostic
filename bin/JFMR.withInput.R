require(remotes)
library(tidyverse) # install.packages("tidyverse")
library(janitor) # install.packages("janitor")
library(scales)
library(mgcv)  # remotes::install_github("ps")

# will ask to install RTools from Rtools35.exe first
library(frs)     # remotes::install_github("ellisp/frs-r-package/pkg", force=T)
library(patchwork)
library(glue) # install.packages("glue") remove.packages("glue")
library(surveillance) # for backprojNP()
library(ggplot2)
library(future)
# Load the package required to read JSON files from web.
library(httr) 
library(jsonlite) 
library(rlist) 
library(dplyr) 
library(data.table)
library(EpiNow2) # remotes::install_github("epiforecasts/EpiNow2")

#### input args
defaultArgs <- list (
    inFile = "WI.rds",
    countyIndex = 1,
    outFile = NULL
)

args <- R.utils::commandArgs(trailingOnly = TRUE,
                             asValues = TRUE ,
                             defaults = defaultArgs)


data_WI_check  <- readRDS(args$inFile)
counties = unique(data_WI_check$NAME)[3:74]

final_df = {}

if (!interactive()){
  options(future.fork.enable = TRUE)
}

#                    *******  CHOOSE YOUR FAVORITE COUNTY HERE  *******
one_county = counties[as.numeric(args$countyIndex)]

data_WI = data_WI_check %>%
    select(NAME, NEGATIVE, POSITIVE, DEATHS, DTH_NEW, POS_NEW, 
           NEG_NEW, TEST_NEW, DATE) %>%
    filter(NAME == one_county) %>% 
    arrange(DATE)%>%
    mutate(date = as.Date(DATE))
  
data_WI = rename(data_WI, confirm = POS_NEW)
data_WI$tests_conducted_total = rowSums(
    cbind(as.numeric(data_WI$NEGATIVE),
          as.numeric(data_WI$POSITIVE)),na.rm=T)
  
d <- data_WI %>%
    replace_na(list(confirm = 0)) %>%
    group_by(date) %>%
    summarise(tests_conducted_total = max(tests_conducted_total, na.rm = TRUE),
              confirm = as.numeric(confirm)) %>%
    mutate(tests_conducted_total  =
               ifelse(tests_conducted_total < 0, NA, tests_conducted_total)
           ) %>%
    ungroup() %>%
    mutate(test_increase = c(
               tests_conducted_total[1], diff(tests_conducted_total)),
           confirm = as.numeric(ifelse(confirm<0, 0, confirm)), 
           pos_raw = pmin(1, confirm / test_increase)) %>%
    complete(date = seq.Date(min(date), max(date), by="day"), 
             fill = list(confirm = 0)) %>%
    mutate(numeric_date = as.numeric(date),
           positivity = pos_raw)

d <- d %>% 
    # in this part the original code was using the fill function for positivity
    # but it's not working for me :( so I filled it with mutate
    mutate(positivity = ifelse( is.nan(positivity) ,
                               mean(d$positivity,na.rm=T), positivity)) %>%
    # I don't believe the sqrt "corrected" cases helped here so have a much more modest 0.1.
    # But first we need to model positivity to smooth it, as it's far too spiky otherwise:
    mutate(ps1 = fitted(gam(positivity ~ s(numeric_date),
                            data = ., family = "quasipoisson")),
           ps2 = fitted(
               loess(positivity ~ numeric_date, data = ., span = 0.1)),
           cases_corrected = confirm * ps1 ^ 0.1 / min(ps1 ^ 0.1)) %>%
    ungroup() %>%
    mutate(smoothed_confirm = fitted(
               loess(confirm ~ numeric_date, data = ., span = 0.1)))

#################################################
# Adjust for testing, delay and right truncation all at once
##################################################
# Sam Abbott, Joel Hellewell, Robin Thompson, Katelyn Gostic, Katharine Sherratt, 
# Sophie Meakin, James Munday, Nikos Bosse and Sebastian Funk (2020). EpiNow2: 
#   Estimate Realtime Case Counts and Time-varying Epidemiological Parameters. 
  
future::plan("multiprocess", gc = TRUE, earlySignal = TRUE)

  
#------------Estimating R with EpiNow2---------------------
# Various delay/lag distributions as per the Covid-19 examples in the EpiNow2 documentation.
  
reporting_delay <- EpiNow2::bootstrapped_dist_fit(rlnorm(100, log(6), 1))
## Set max allowed delay to 30 days to truncate computation
reporting_delay$max <- 30

generation_time <- list(mean = EpiNow2::covid_generation_times[1, ]$mean,
                        mean_sd = EpiNow2::covid_generation_times[1, ]$mean_sd,
                        sd = EpiNow2::covid_generation_times[1, ]$sd,
                        sd_sd = EpiNow2::covid_generation_times[1, ]$sd_sd,
                        max = 30)

incubation_period <- list(mean = EpiNow2::covid_incubation_period[1, ]$mean,
                          mean_sd = EpiNow2::covid_incubation_period[1, ]$mean_sd,
                          sd = EpiNow2::covid_incubation_period[1, ]$sd,
                          sd_sd = EpiNow2::covid_incubation_period[1, ]$sd_sd,
                          max = 30)

#future::plan(sequential)
#---------Based on positivity-adjusted-------------

d2 <- select(d, date, cases_corrected) %>%
  mutate(confirm = round(cases_corrected) )

#future::plan("multiprocess", gc = TRUE, earlySignal = TRUE)

estimates2 <- EpiNow2::epinow(reported_cases = d2, 
                              generation_time = generation_time,
                              delays = list(incubation_period, reporting_delay),
                              horizon = 7, samples = 35, warmup = 10, 
                              cores = 1, chains = 1, verbose = TRUE, 
                              adapt_delta = 0.95)

future::plan(sequential)

estimation_data = estimates2$plots$summary$data
estimation_data$county = rep(one_county, dim(estimation_data)[1] )

final_df = rbind(final_df, estimation_data)

write.table(final_df, args$outFile,
            col.names=F,row.names = FALSE, 
            sep = "\t", quote = FALSE)

q()
