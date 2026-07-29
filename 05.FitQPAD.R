library(tidyverse)
library(RTMB)
library(terra)
library(sf)
library(foreach)
library(doParallel)

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
root <- getwd() # if not running on cluster, comment this out

v.wt <- "2026-07-10"
load(file.path(root, "WildTrax", v.wt, paste0("02_wildtrax_clean_", v.wt, ".Rdata")))

qpad_dir_out = file.path(root, "WildTrax", v.wt, "qpad_fits")
if (!dir.exists(qpad_dir_out)) dir.create(qpad_dir_out)

MAX_DIST = 1000 # in m
R_SCALE = 1 / 1000 # m to km
T_SCALE = 1 / 60 # seconds to minutes
DEG_CRS = "EPSG:4326"
INITS = c(0, 1)
MIN_TOD_OBS = 5 # minimum # of observations for a given time of day for the models to work; if there are fewer than 5 observations for a given time of day category, those will be removed
DEG_DIGITS = 4 # how to round degrees, I know it's not an exact spatial unit but close enough; this is equal to about 100 m

all_species = unique(pc.good.final$species_code)
if (!is.na(commandArgs()[6])) all_species = str_split_1(commandArgs()[6], " ")

# message("length(commandArgs()) = ", length(commandArgs()))
# print(commandArgs())

lambda_covs_formula = ~ 0 + morning + sunrise + day + sunset + night + nauticaldawn + nauticaldusk + is_pointcount + timeofyear
alpha_covs_formula = ~ open_closed

# Prepare PC and ARU data for RTMB models - do it now so it's a bit faster ----
aru_final = aru.good %>%
  dplyr::filter(!is.na(individual_count)) %>%
  mutate(r_lo = 0,
         r_up = MAX_DIST * R_SCALE,
         r_max = MAX_DIST * R_SCALE,
         t_lo = detection_time * T_SCALE,
         t_up = detection_time * T_SCALE,
         t_max = task_duration * T_SCALE,
         count = individual_count,
         r_mxs = r_max,
         type = "aru") %>%
  dplyr::select(project_id, survey_id = task_id, species_code, r_lo, r_up, r_max, t_lo, t_up, t_max, count, r_mxs, longitude, latitude, date = recording_date_time, type)

pc_final = pc.good.final %>%
  dplyr::filter(!is.na(individual_count)) %>%
  # correct names for distance and duration methods
  mutate(detection_distance = case_match(detection_distance,
                                         "UNKNOWN" ~ "0m-INF",
                                         "0m-INF_ARU" ~ "0m-INF",
                                         NA ~ "0m-INF",
                                         .default = detection_distance),
         survey_distance_method = ifelse(detection_distance == "0m-INF", "0m-INF", survey_distance_method),
         survey_distance_method = case_match(survey_distance_method,
                                             "UNKNOWN" ~ "0m-INF",
                                             "0m-INF_ARU" ~ "0m-INF",
                                             NA ~ "0m-INF",
                                             .default = survey_distance_method),
         # replace all "INF" with the maximum distance 
         survey_distance_method = str_replace_all(survey_distance_method, "INF", paste0(MAX_DIST, "m")),
         detection_distance = str_replace_all(detection_distance, "INF", paste0(MAX_DIST, "m")),
         survey_duration_method = ifelse(detection_time == "UNKNOWN", paste0("0-", str_split_i(survey_duration_method, "-", -1)), survey_duration_method),
         detection_time = ifelse(detection_time == "UNKNOWN", paste0("0-", str_split_i(survey_duration_method, "-", -1)), detection_time),
         survey_duration_method = case_match(survey_duration_method,
                                             "0-3-5-10-10min+" ~ "0-3-5-10min",
                                             .default = survey_duration_method),
         n_distance_bins = str_count(survey_distance_method, "-"),
         n_duration_bins = str_count(survey_duration_method, "-"),
         n_distance_bins = ifelse(is.na(n_distance_bins), 1, n_distance_bins),
         n_duration_bins = ifelse(is.na(n_duration_bins), 1, n_duration_bins),
         r_max = as.numeric(str_split_i(survey_distance_method, "(-|[a-z])+", -2)),
         t_max = as.numeric(str_split_i(survey_duration_method, "(-|[a-z])+", -2)),
         r_lo = as.numeric(str_split_i(detection_distance, "(-|[a-z])+", 1)),
         r_up = as.numeric(str_split_i(detection_distance, "(-|[a-z])+", -2)),
         t_lo = as.numeric(str_split_i(detection_time, "(-|[a-z])+", 1)),
         t_up = as.numeric(str_split_i(detection_time, "(-|[a-z])+", -2))) %>% 
  # remove sightings without any duration bins
  dplyr::filter(n_distance_bins > 1 | n_duration_bins > 1) %>%
  dplyr::select(project_id, survey_id, species_code, r_lo, r_up, r_max, t_lo, t_up, t_max, count = individual_count, longitude, latitude, date = survey_date) %>%
  mutate(r_lo = r_lo * R_SCALE,
         r_up = r_up * R_SCALE,
         r_mxs = r_max * R_SCALE,
         r_max = MAX_DIST * R_SCALE,
         type = "pc") 

all_final = rbind(pc_final, aru_final) %>%
  group_by(species_code) %>%
  mutate(n_this_sp = n()) %>%
  ungroup %>%
  arrange(-n_this_sp, species_code)

rm(pc.good.final, aru.good, pc_final, aru_final)
gc()

# Extract covariates ----
timeofday_cov = cov_tod_bin("timeofday")

closed_cov = cov_dynamic_raster(cov_name = "open_closed", raster_dir = file.path("data", "scanfi_biomass_agg"), method = "nearest")

all_final = all_final %>%
  mutate(timeofday = timeofday_cov$get(longitude, latitude, date, crs_in = DEG_CRS),
         open_closed = closed_cov$get(longitude, latitude, date, crs_in = DEG_CRS))

wt.wide = wt.wide %>%
  mutate(timeofday = timeofday_cov$get(longitude, latitude, date, crs_in = DEG_CRS),
         open_closed = closed_cov$get(longitude, latitude, date, crs_in = DEG_CRS))

# Set up parallelization
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) * as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
if (is.na(ncores)) ncores = 5
registerDoParallel(cores = ncores)

qpad_fits = foreach(sp = all_species, .errorhandling = "stop") %dopar% {
  
  qpad_fit_dir_out = file.path(qpad_dir_out, paste0(sp, ".rds"))
  
  if (file.exists(qpad_fit_dir_out)) {
    message(sp, " is already done. ", Sys.time())
    return(0)
  }
  
  message("beginning ", sp, ": ", Sys.time())
  
  this_detections = all_final %>%
    dplyr::filter(species_code == sp) %>%
    group_by(survey_id) %>%
    mutate(count_total = sum(count)) %>%
    ungroup %>%
    mutate(lon_rounded = round(longitude, DEG_DIGITS),
           lat_rounded = round(latitude, DEG_DIGITS),
           year = year(date)) %>%
    group_by(lon_rounded, lat_rounded, year) %>%
    mutate(max_count = max(count_total)) %>%
    ungroup
  
  # Find the surveys at each location-year combo that returned the max count for the species, and then figure out how long those surveys lasted
  rounded_duration_info = this_detections %>%
    group_by(lon_rounded, lat_rounded, year) %>%
    dplyr::filter(count_total == max_count) %>%
    summarize(t_mxm = min(t_max),
              best_survey_id = survey_id[which.min(t_max)])
  
  this_detections = this_detections %>%
    left_join(rounded_duration_info, by = join_by(lon_rounded == lon_rounded, lat_rounded == lat_rounded, year == year))
  
  rounded_survey_info = this_detections %>%
    group_by(lon_rounded, lat_rounded, year) %>%
    summarize(n_surveys = n_distinct(survey_id),
              max_count = max(max_count)) %>%
    left_join(rounded_duration_info, by = join_by(lon_rounded == lon_rounded, lat_rounded == lat_rounded, year == year))
  
  # Add surveys with 0 birds observed for additional right censored samples
  wt.wide$this_species = wt.wide[, sp, drop = TRUE] # so we have a consistent column name for upcoming dplyr operations
  
  rc_surveys = wt.wide %>%
    mutate(lon_rounded = round(longitude, DEG_DIGITS),
           lat_rounded = round(latitude, DEG_DIGITS),
           year = year(date_time),
           .after = distance) %>%
    left_join(rounded_survey_info, by = join_by(lon_rounded == lon_rounded, lat_rounded == lat_rounded, year == year)) %>%
    dplyr::filter(n_surveys > 0, this_species == 0)
  
  # align formatting with all_rtmb_ready
  all_rtmb_rc = rc_surveys %>%
    mutate(distance = ifelse(is.finite(distance), distance, MAX_DIST) / 1000,
           duration = duration / 60,
           obs_id = paste0("RC_", seq_len(n())),
           r_lo = 0,
           r_up = distance,
           r_max = distance,
           t_lo = duration,
           t_up = Inf,
           t_max = duration,
           count = 0,
           count_total = 0,
           r_mxs = distance,
           type = ifelse(method == "PC", "pc", "aru"),
           equipment_model = NA) %>%
    dplyr::select(obs_id, survey_id, best_survey_id, project_id, r_lo, r_up, r_max, t_lo, t_up, t_max, t_mxm, count, r_mxs, longitude, latitude, date = date_time, type, task_method = method, equipment_model, count_total, lon_rounded, lat_rounded, year, max_count)
  
  all_rtmb_ready = rbind(this_detections, all_rtmb_rc) %>%
    group_by(timeofday) %>%
    mutate(n_this_tod = sum(pmin(count, 1))) %>%
    ungroup %>%
    dplyr::filter(n_this_tod >= MIN_TOD_OBS) %>%
    mutate(nauticaldawn = timeofday == "nauticaldawn",
           sunrise = timeofday == "sunrise",
           morning = timeofday == "morning",
           day = timeofday == "day",
           sunset = timeofday == "sunset",
           nauticaldusk = timeofday == "nauticaldusk",
           night = timeofday == "night",
           open_closed = open_closed / max(open_closed),
           timeofyear = yday(date) - 152,
           is_pointcount = (type == "pc"))
  
  obj = try(fit_jqpadrc(all_rtmb_ready, 
                        formula_alpha = alpha_covs_formula, 
                        formula_lambda = lambda_covs_formula,
                        return_data = TRUE, 
                        inits = INITS, 
                        profile_improve_stop = 1, 
                        return_hess = TRUE, 
                        return_ci = TRUE))

  message("completed ", sp, ": ", Sys.time())
  saveRDS(list(full = obj), qpad_fit_dir_out)
  
  obj
  
}

error_fits = !sapply(qpad_fits, is.numeric)
if (any(error_fits)) message("Fits that returned errors: [", paste0(all_species[error_fits], collapse = ", "), "]")

# names(qpad_fits) = all_species
# save(qpad_fits, file = file.path(root, "WildTrax", v.wt, paste0("05_qpad_estimates_", v.wt, ".Rdata")))
