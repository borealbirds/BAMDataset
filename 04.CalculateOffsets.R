library(tidyverse)
library(RTMB)
library(sf)
library(duckdb)
library(duckspatial)

# Read in BAM Dataset ----

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
v.wt <- "2026-08-26"
v.ebd = "Jun-2026"
v.prism = "2026-10-09"

# Read in QPAD data and filter to "useful" observations with distance / time info ----

db_conn_path = file.path(root, paste0("03_QPADDataset_WT-", v.wt, "_EBd-", v.ebd, "_PRISM-", v.prism, ".duckdb"))
db_conn = dbConnect(duckdb(db_conn_path))
dbExecute(db_conn, "LOAD spatial") # have to do this for whatever reason

db_bam_conn_path = file.path(root, paste0("03_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, "_PRISM-", v.prism, ".duckdb"))
db_bam_conn = dbConnect(duckdb(db_bam_conn_path))
dbExecute(db_bam_conn, "LOAD spatial") # have to do this for whatever reason

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