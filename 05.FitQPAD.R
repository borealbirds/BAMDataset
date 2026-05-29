library(tidyverse)
library(RTMB)
library(foreach)
library(doParallel)

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

root <- "inputs"
# root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

v.wt <- "2026-05-28"
load(file.path(root, "WildTrax", v.wt, paste0("02_wildtrax_clean_", v.wt, ".Rdata")))

MAX_DIST = 1000 # in m
R_SCALE = 1 / 1000 # m to km
T_SCALE = 1 / 60 # seconds to minutes

all_species = unique(pc.good.final$species_code)

# Set up parallelization
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) * as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
if (is.na(ncores)) ncores = 5
registerDoParallel(cores = ncores)

qpad_fits = foreach(sp = all_species) %dopar% {
  
  aru_this = aru.good %>% dplyr::filter(species_code == SP, !is.na(individual_count))
  pc_this = pc.good.final %>% dplyr::filter(species_code == SP, !is.na(individual_count))
  
  pc_rtmb_ready = pc_this %>%
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
    dplyr::select(r_lo, r_up, r_max, t_lo, t_up, t_max, count = individual_count, longitude, latitude, date = survey_date) %>%
    mutate(r_lo = r_lo * R_SCALE,
           r_up = r_up * R_SCALE,
           r_mxs = r_max * R_SCALE,
           r_max = MAX_DIST * R_SCALE,
           type = "pc") 
  
  aru_rtmb_ready = aru_this %>%
    mutate(r_lo = 0,
           r_up = MAX_DIST * R_SCALE,
           r_max = MAX_DIST * R_SCALE,
           t_lo = detection_time * T_SCALE,
           t_up = detection_time * T_SCALE,
           t_max = task_duration * T_SCALE,
           count = individual_count,
           r_mxs = r_max,
           type = "aru") %>%
    dplyr::select(r_lo, r_up, r_max, t_lo, t_up, t_max, count, r_mxs, longitude, latitude, date = recording_date_time, type)
  
  all_rtmb_ready = rbind(pc_rtmb_ready, aru_rtmb_ready)
  
  # Extract covariates
  closed_cov = cov_dynamic_raster(cov_name = "open_closed",
                                  raster_dir = file.path("data", "scanfi_biomass_agg"),
                                  method = "nearest")
  t_since_rise_cov = cov_timeofday(cov_name = "t_since_sunrise",
                                   type = "sunrise")
  t_since_set_cov = cov_timeofday(cov_name = "t_since_sunset",
                                  type = "sunset")
  
  closed_cov_values = with(all_rtmb_ready, closed_cov$get(longitude, latitude, date, crs_in = "EPSG:4326"))
  sunrise_values = with(all_rtmb_ready, t_since_rise_cov$get(longitude, latitude, date, crs_in = "EPSG:4326"))
  sunset_values = with(all_rtmb_ready, t_since_set_cov$get(longitude, latitude, date, crs_in = "EPSG:4326"))
  
  all_rtmb_ready = all_rtmb_ready %>% 
    mutate(open_closed = closed_cov_values / max(closed_cov_values),
           t_since_sunrise = sunrise_values,
           t_since_sunset = sunset_values,
           dawn = t_since_sunrise >= -0.5 & t_since_sunrise < 2,
           midday = t_since_sunrise >= 2 & t_since_sunset < 0.5,
           dusk = t_since_sunset >= -0.5 & t_since_sunset < 2,
           night = !(dawn | midday | dusk))
  
  pc_rtmb_ready_cov = all_rtmb_ready %>% dplyr::filter(type == "pc")
  aru_rtmb_ready_cov = all_rtmb_ready %>% dplyr::filter(type == "aru")
  
  lambda_covs_formula = ~ 0 + night + dawn + midday + dusk
  alpha_covs_formula = ~ open_closed
  
  obj_pc_null = try(fit_jqpadmix(pc_rtmb_ready_cov, return_data = TRUE))
  # obj_aru_null = try(fit_jqpadmix(aru_rtmb_ready_cov, return_data = TRUE))
  # obj_null = fit_jqpadmix(all_rtmb_ready, return_data = TRUE)
  
  obj_pc = try(fit_jqpadmix(pc_rtmb_ready_cov, formula_alpha = alpha_covs_formula, formula_lambda = lambda_covs_formula, return_data = TRUE))
  # obj_aru = try(fit_jqpadmix(aru_rtmb_ready_cov, formula_alpha = alpha_covs_formula, formula_lambda = lambda_covs_formula, return_data = TRUE))
  # obj = fit_jqpadmix(all_rtmb_ready, formula_alpha = alpha_covs_formula, formula_lambda = lambda_covs_formula, return_data = TRUE)
  
  list(full = obj_pc, null = obj_pc_null)
  
}

names(qpad_fits) = all_species

save(qpad_fits, file = file.path(root, "WildTrax", v.wt, paste0("05_qpad_estimates_", v.wt, ".Rdata")))
