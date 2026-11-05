library(tidyverse)
library(RTMB)
library(sf)
library(duckdb)
library(duckspatial)
library(terra)
library(mvtnorm)
library(dtplyr)

# Read in BAM Dataset ----

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

FIR = str_detect(getwd(), "^/project")
if (!FIR) {
  root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
  qpad_report_dir = "G:/Shared drives/BAM_RshProjs/Methodology/Detectability/QPADModelFits/qpad_fits"
  qpad_fit_dir = file.path(qpad_report_dir, c("Other", "Typicaldawnspecies", "Waterbirds"))
  dat_dir = "G:/Shared drives/BAM_RshProjs/Methodology/Detectability/QPADModelFits/data"
} else {
  root = getwd()
  qpad_fit_dir = "/scratch/pt1/bam/qpad_fits"
  qpad_report_dir = qpad_fit_dir
  dat_dir = "data"
}

search_for_file = function(dirs, filename) {
  
  all_files = file.path(dirs, filename)
  all_files[file.exists(all_files)]
  
}

v.wt <- "2026-08-26"
v.ebd = "Jun-2026"
v.prism = "2026-10-09"

duckdb_lock_file = "duckdb_bam_lock.csv"

# if there's already a lock file, just run this again and again and wait. Add an extra Sys.sleep to space out in case a bunch of jobs make the file at the same time
Sys.sleep(runif(1, 0, 30))
while(file.exists(duckdb_lock_file)) Sys.sleep(0.5)

# once the file no longer exists, write it
write_csv(data.frame(0), duckdb_lock_file)

# Read in QPAD data and filter to "useful" observations with distance / time info ----

db_conn_path = file.path(root, paste0("03_QPADDataset_WT-", v.wt, "_EBd-", v.ebd, "_PRISM-", v.prism, ".duckdb"))
db_conn = dbConnect(duckdb(db_conn_path))
dbExecute(db_conn, "LOAD spatial") # have to do this for whatever reason

db_bam_conn_path = file.path(root, paste0("03_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, "_PRISM-", v.prism, ".duckdb"))
db_bam_conn = dbConnect(duckdb(db_bam_conn_path))
dbExecute(db_bam_conn, "LOAD spatial") # have to do this for whatever reason

visit = ddbs_read_table(db_bam_conn, name = "visit")

# disconnect!
dbDisconnect(db_conn)
dbDisconnect(db_bam_conn)

# remove the lock so other processes can run
file.remove(duckdb_lock_file)

# Add covariates ----

DEG_CRS = "EPSG:4326"
N_BOOT = 30 # for generating random offsets

timeofday_cov = cov_tod_bin("timeofday")
closed_cov = cov_dynamic_raster(cov_name = "open_closed", raster_dir = file.path(dat_dir, "scanfi_biomass_agg"), method = "nearest")
road_cov = cov_static_raster(cov_name = "near_road", raster_path = file.path(dat_dir, "all_roads_buffered.tif"), na_in_domain = TRUE)

visit_cov = visit[1:1000, ] %>%
  dplyr::filter(!is.na(longitude), !is.na(latitude), !is.na(date_time)) %>%
  mutate(timeofday = timeofday_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS),
         biomass = closed_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS) / 200,
         near_road = road_cov$get(longitude, latitude, date_time, crs_in = DEG_CRS),
         timeofyear = yday(date_time) - 152,
         aru_01 = (sensor == "ARU"),
         distance_km = pmin(1, distance / 1000),
         duration_min = duration / 60) %>%
  dplyr::filter(!is.na(timeofday), !is.na(distance_km), !is.na(duration_min))

pars_df_base = expand_grid(boot_id = seq_len(N_BOOT),
                           survey_id = visit_cov$survey_id) %>%
  left_join(visit_cov, by = join_by(survey_id == survey_id))

# Go through each species and calculate offset ----

ALL_POSSIBLE_PAR_NAMES = c("log_lambda_nauticaldawn",
                           "log_lambda_sunrise",
                           "log_lambda_morning",
                           "log_lambda_day",
                           "log_lambda_sunset",
                           "log_lambda_nauticaldusk",
                           "log_lambda_night",
                           "log_lambda_timeofyear",
                           "log_lambda_is_aru",
                           "log_alpha_",
                           "log_alpha_open_closed",
                           "log_alpha_near_road")

all_species_df = read_csv(file.path(qpad_report_dir, "reports", "qpad_report_20261028.csv"))
all_species = all_species_df %>%
  dplyr::filter(!is.na(species_code)) %>%
  pull(species_code)

all_offsets = lapply(all_species[1:3], function(sp) {
  
  message("beginning ", sp)
  
  this_report = all_species_df %>%
    dplyr::filter(species_code == sp)
  this_report_timeofday = this_report %>%
    dplyr::select(include_day:include_sunset) %>%
    pivot_longer(cols = everything(), names_to = "timeofday", values_to = "include", names_prefix = "include_")
  tod_to_keep = this_report_timeofday %>%
    dplyr::filter(include) %>%
    pull(timeofday)
  include_this_species = any(this_report_timeofday$include)
  include_road = this_report$include_road
  include_aru = this_report$include_aru
  
  if (!is.na(include_this_species) && include_this_species) {
    ## Read in the appropriate QPAD fit ----
    
    this_fit_filepath = search_for_file(qpad_fit_dir, paste0(sp, ".rds"))
    if (!file.exists(this_fit_filepath)) {
      message("There is no existing QPAD fit for ", sp, ", even though it is in the report. Using an old report file? Moving on...")
    } else {
      
      message("Reading in the fit...")
      
      this_fit = readRDS(this_fit_filepath)$full
      
      ## Calculate the offset for each survey ----
      this_pars_all = this_fit$fit$par
      lambda_par_names = paste0("log_lambda_", colnames(this_fit$dat$lambda_covs))
      alpha_par_names = paste0("log_alpha_", colnames(this_fit$dat$alpha_covs))
      
      message("Making random parameters...")
      
      r_pars = rmvnorm(N_BOOT, this_pars_all, get_cov_matrix(this_fit)) %>%
        as.data.frame %>%
        set_names(c(lambda_par_names, alpha_par_names))
      
      pars_df = cbind(pars_df_base, r_pars[rep(seq_along(nrow(r_pars)), each = N_BOOT), ]) %>%
        check_names(ALL_POSSIBLE_PAR_NAMES, fill = 0) %>%
        mutate(log_lambda = log_lambda_nauticaldawn * (timeofday == "nauticaldawn") +
                 log_lambda_sunrise * (timeofday == "sunrise") +
                 log_lambda_morning * (timeofday == "morning") +
                 log_lambda_day * (timeofday == "day") +
                 log_lambda_sunset * (timeofday == "sunset") +
                 log_lambda_nauticaldusk * (timeofday == "nauticaldusk") +
                 log_lambda_night * (timeofday == "night") +
                 log_lambda_timeofyear * timeofyear,
               log_alpha = log_alpha_ + 
                 log_alpha_near_road * ifelse(!is.na(include_road) && include_road, near_road, 0) + 
                 log_alpha_open_closed * biomass)
      
      prism_conversion_rates = log(pi * visit_cov$distance_km ^ 2 * 100) # area of PRISM plots in hectares
      
      message("Calculating offsets and metrics...")
      
      all_par_offsets = pars_df %>%
        mutate(offset = get_offset_la(log_lambda, log_alpha, duration_min, distance_km),
               # fix if NA
               offset = ifelse(timeofday %in% tod_to_keep, offset, NA)) %>%
        group_by(survey_id) %>%
        # placeholder names for now so we can use dplyr notation
        summarize(mean_log_offset = mean(offset),
                  sd_log_offset = sd(offset)) %>%
        ungroup %>%
        # correct for PRISM
        mutate(sensor = visit_cov$sensor,
               mean_log_offset = ifelse(sensor == "PRISM", prism_conversion_rates, mean_log_offset),
               sd_log_offset = ifelse(sensor == "PRISM", 0, sd_log_offset)) %>%
        # remove columns that would otherwise be repeated when we cbind
        dplyr::select(-survey_id, -sensor) %>%
        # species-specific names to use in the final join
        set_names(paste0("offset_", c("mean", "sd"), "_", sp))
      
      message("Returning...")
      
      return(all_par_offsets)
    }
  }
  
  numeric(0)
  
})

all_offsets_df = do.call(cbind, all_offsets[sapply(all_offsets, is.data.frame)]) %>%
  mutate(survey_id = visit_cov$survey_id)

# Add new dataframe with offsets to duckdb connection and write ----

# if there's already a lock file, just run this again and again and wait. I don't think we need the extra Sys.sleep() anymore.
while(file.exists(duckdb_lock_file)) Sys.sleep(0.5)

# once the file no longer exists, write it
write_csv(data.frame(0), duckdb_lock_file)

db_bam_conn = dbConnect(duckdb(db_bam_conn_path))
dbExecute(db_bam_conn, "LOAD spatial") # have to do this for whatever reason

dbWriteTable(db_bam_conn, name = "offsets", all_offsets_df, overwrite = TRUE)
dbDisconnect(db_bam_conn)

# remove the lock so other processes can run
file.remove(duckdb_lock_file)
