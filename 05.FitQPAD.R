library(tidyverse)
library(RTMB)
library(terra)
library(sf)
library(foreach)
library(doParallel)

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
root <- getwd() # if not running on cluster, comment this out

v.wt <- "2026-06-02"
load(file.path(root, "WildTrax", v.wt, paste0("02_wildtrax_clean_", v.wt, ".Rdata")))

qpad_dir_out = file.path(root, "WildTrax", v.wt, "qpad_fits")
if (!dir.exists(qpad_dir_out)) dir.create(qpad_dir_out)

MAX_DIST = 1000 # in m
R_SCALE = 1 / 1000 # m to km
T_SCALE = 1 / 60 # seconds to minutes
DEG_CRS = "EPSG:4326"
INITS = c(0, 1)

all_species = unique(pc.good.final$species_code)
# all_species = c("WOTH", "GCTH")

lambda_covs_formula = ~ 0 + dawn + morning + midday + dusk + night + highlat
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
  dplyr::select(species_code, r_lo, r_up, r_max, t_lo, t_up, t_max, count, r_mxs, longitude, latitude, date = recording_date_time, type)

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
  dplyr::select(species_code, r_lo, r_up, r_max, t_lo, t_up, t_max, count = individual_count, longitude, latitude, date = survey_date) %>%
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

rm(pc.good.final, wt.wide, aru.good, pc_final, aru_final)
gc()

# Extract covariates ----
closed_cov = cov_dynamic_raster(cov_name = "open_closed", raster_dir = file.path("data", "scanfi_biomass_agg"), method = "nearest")
t_since_rise_cov = cov_timeofday(cov_name = "t_since_sunrise", type = "sunrise")
t_since_set_cov = cov_timeofday(cov_name = "t_since_sunset", type = "sunset")

all_final$open_closed = with(all_final, closed_cov$get(longitude, latitude, date, crs_in = DEG_CRS))
all_final$t_since_sunrise = with(all_final, t_since_rise_cov$get(longitude, latitude, date, crs_in = DEG_CRS))
all_final$t_since_sunset = with(all_final, t_since_set_cov$get(longitude, latitude, date, crs_in = DEG_CRS))

# Set up parallelization
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) * as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
if (is.na(ncores)) ncores = 5
registerDoParallel(cores = ncores)

qpad_fits = foreach(sp = all_species, .errorhandling = "stop") %dopar% {
  
  qpad_fit_dir_out = file.path(qpad_dir_out, paste0(sp, ".rds"))
  
  if (file.exists(qpad_fit_dir_out)) {
    message(sp, " is already done. ", Sys.time())
  } else {
    message("beginning ", sp, ": ", Sys.time())
    
    all_rtmb_ready = all_final %>%
      dplyr::filter(species_code == sp) %>% 
      mutate(t_since_sunrise = ifelse(is.na(t_since_sunrise), Inf, t_since_sunrise),
             t_since_sunset = ifelse(is.na(t_since_sunset), Inf, t_since_sunset),
             dawn = t_since_sunrise >= -0.5 & t_since_sunrise < 1,
             morning = t_since_sunrise >= 1 & t_since_sunrise < 5,
             midday = t_since_sunrise >= 5 & t_since_sunset < -1,
             dusk = t_since_sunset >= -1 & t_since_sunset < 0.5,
             night = !(dawn | midday | morning | dusk),
             highlat = is.infinite(t_since_sunrise) | is.infinite(t_since_sunset))
    
    obj_null = fit_jqpadmix(all_rtmb_ready, return_data = TRUE, inits = INITS, profile_improve_stop = 1)
    obj = fit_jqpadmix(all_rtmb_ready, formula_alpha = alpha_covs_formula, formula_lambda = lambda_covs_formula, return_data = TRUE, inits = INITS, profile_improve_stop = 1)
    
    message("completed ", sp, ": ", Sys.time())
    saveRDS(list(full = obj, null = obj_null), qpad_fit_dir_out)
  }
  0
  # qpad_fits = list(full = obj_pc, null = obj_pc_null)
  
}

error_fits = !sapply(qpad_fits, is.numeric)
if (any(error_fits)) message("Fits that returned errors: [", paste0(all_species[error_fits], collapse = ", "), "]")

# names(qpad_fits) = all_species
# save(qpad_fits, file = file.path(root, "WildTrax", v.wt, paste0("05_qpad_estimates_", v.wt, ".Rdata")))
