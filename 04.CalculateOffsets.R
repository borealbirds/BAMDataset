library(tidyverse)
library(RTMB)
library(sf)
library(duckdb)
library(duckspatial)
library(terra)

# Read in BAM Dataset ----

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

FIR = str_detect(getwd(), "^/project")
if (!FIR) {
  root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
  qpad_fit_dir = "G:/Shared drives/BAM_RshProjs/Methodology/Detectability/QPADModelFits/qpad_fits"
  dat_dir = "G:/Shared drives/BAM_RshProjs/Methodology/Detectability/QPADModelFits/data"
} else {
  root = ""
  qpad_fit_dir = "/scratch/pt1/bam/qpad_fits"
  dat_dir = "data"
}

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

visit = ddbs_read_table(db_bam_conn, name = "visit")

# Add covariates ----

DEG_CRS = "EPSG:4326"

timeofday_cov = cov_tod_bin("timeofday")
closed_cov = cov_dynamic_raster(cov_name = "open_closed", raster_dir = file.path(dat_dir, "scanfi_biomass_agg"), method = "nearest")
road_cov = cov_static_raster(cov_name = "near_road", raster_path = file.path(dat_dir, "all_roads_buffered.tif"), na_in_domain = TRUE)

visit_cov = visit %>%
  dplyr::filter(!is.na(longitude), !is.na(latitude), !is.na(date_time)) %>%
  mutate(timeofday = timeofday_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS),
         biomass = closed_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS) / 200,
         near_road = road_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS),
         timeofyear = yday(date_time) - 152,
         aru_01 = (sensor == "ARU"),
         distance_km = pmin(1, distance / 1000),
         duration_min = duration / 60) %>%
  dplyr::filter(!is.na(timeofday))

# Go through each species and calculate offset ----

all_species_df = read_csv(file.path(qpad_fit_dir, "qpad_report_20261013.csv"))
all_species = all_species_df %>%
  dplyr::filter(!is.na(species_code)) %>%
  pull(species_code)

for (sp in all_species) {
  
  this_report = all_species_df %>%
    dplyr::filter(species_code == sp)
  this_report_timeofday = this_report %>%
    dplyr::select(include_day:include_sunset) %>%
    pivot_longer(cols = everything(), names_to = "timeofday", values_to = "include", names_prefix = "include_")
  include_this_species = any(this_report_timeofday$include)
  include_road = this_report$include_road
  include_aru = this_report$include_aru
  
  if (!is.na(include_this_species) && include_this_species) {
    ## Read in the appropriate QPAD fit ----
    
    this_fit_filepath = file.path(qpad_fit_dir, paste0(sp, ".rds"))
    this_fit = readRDS(this_fit_filepath)$full
    
    ## Calculate the offset for each survey ----
    this_pars_all = this_fit$fit$par
    this_lambda_pars = this_pars_all[str_detect(names(this_pars_all), "^log_lambda_")]
    this_lambda_pars_named = this_lambda_pars %>%
      set_names(paste0("log_lambda_", colnames(this_fit$dat$lambda_covs)))
    this_alpha_pars = this_pars_all[str_detect(names(this_pars_all), "^log_alpha_")]
    this_alpha_pars_named = this_alpha_pars %>%
      set_names(paste0("log_alpha_", colnames(this_fit$dat$alpha_covs)))
    
    these_offsets = with(visit_cov, get_offset(this_lambda_pars_named,
                                               this_alpha_pars_named,
                                               duration_min,
                                               distance_km,
                                               timeofday,
                                               timeofyear,
                                               biomass,
                                               ifelse(!is.na(include_road) && include_road, near_road, 0),
                                               # assume all counts are point counts for purposes of offsets
                                               0))
    
    prism_conversion_rates = pi * visit_cov$distance_km ^ 2 * 100 # area of PRISM plots in hectares
    new_offset_df = data.frame(offset = ifelse(visit_cov$sensor == "PRISM", prism_conversion_rates, these_offsets)) %>%
      cbind(data.frame(timeofday = visit_cov$timeofday)) %>%
      left_join(this_report_timeofday, by = join_by(timeofday == timeofday)) %>%
      mutate(offset = ifelse(include, offset, NA)) %>%
      pull(offset) %>%
      as.data.frame %>%
      set_names(paste0("offset_", sp))
    
    visit_cov = cbind(visit_cov, new_offset_df)
  }
  
}

# Add new dataframe with offsets to duckdb connection and write ----

offsets = visit_cov %>%
  st_drop_geometry %>%
  dplyr::select(survey_id, starts_with("offset_"))

dbWriteTable(db_bam_conn, name = "offsets", offsets)
