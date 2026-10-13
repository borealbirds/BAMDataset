library(tidyverse)
library(RTMB)

# Read in BAM Dataset ----



# Add covariates ----



# Go through each species and calculate offset ----

get_offset = function(lambda_pars,
                      alpha_pars,
                      duration = 3,
                      distance = 1,
                      timeofday = "sunrise",
                      timeofyear = 0,
                      biomass = 0,
                      near_road = 0,
                      aru_01 = 0) {
  
  
  
}

for (sp in all_species) {
  
  ## Read in the appropriate QPAD fit ----
  
  this_fit_filepath = file.path(qpad_fit_dir, paste0(sp, ".rds"))
  this_fit = readRDS(this_fit_filepath)
  
  ## Calculate the offset for each survey ----
  
  these_offsets = with(visit, get_offset(this_lambda_pars,
                                         this_alpha_pars,
                                         duration,
                                         distance,
                                         timeofday,
                                         timeofyear,
                                         biomass,
                                         near_road,
                                         aru_01))
  
  visit[[paste0("offset_", sp)]] = these_offsets
  
}

# Add new dataframe with offsets to duckdb connection and write ----