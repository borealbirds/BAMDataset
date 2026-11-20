# ---
# title: BAM dataset - harmonize data
# author: Elly Knight
# created: June 10, 2026
# ---

#NOTES################################

#PURPOSE: This script removes datasets that are only to be used in national models to produce a dataset for other BAM objectives. These are datasets identified by partners for removal unless prior approval is received. Currently these projects only include ones from CWS-NOR.

#TODO: Change to 05_ object & updated when it's ready

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Set the WildTrax version ----
v.wt <- "2026-06-02"

#4. Set the eBird version ----
v.ebd <- "Jan-2026"

#REMOVE DATASETS#########

#1. Get dataset ----
load(file.path(root, paste0("04_BAMDataset_WT-", v.wt, "_EBd-", v.ebd,  ".Rdata")))

#2. Get project details from WildTrax download ----
load(file.path(root, "WildTrax", v.wt, paste0("01_wildtrax_raw_", v.wt, ".Rdata")))

#3. Identify project ids to remove ----
# As per communication from CWS-NOR (Eamon Riordan-Short) and CWS-ON (Kevin Hannah)
ids <- proj |> 
  dplyr::filter(organization_name=="CWS-NOR" |
                organization_name=="CWS-ONT" & str_detect(project, "HADMU")) |> 
  dplyr::filter(project_status!="Published - Public")

#4. Filter ----
dat <- dat |> 
  dplyr::filter(!project_id %in% ids$project_id)

#5. Save ----
save(dat, file=file.path(root, paste0("07_BAMDataset_NonModels_WT-", v.wt, "_EBd-", v.ebd,  ".Rdata")))
