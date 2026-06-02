# ---
# title: BAM dataset - clean WildTrax data
# author: Elly Knight
# created: March 2, 2026
# ---

#NOTES################################

#PURPOSE: This script tidies the WT data downloads into objects that can then be harmonized with eBird data to produce the final dataset.

#Note that surveys done with acoustic recorders and transcribed to point count and stored in the point count sensor (survey_distance_method=="0m-INF-ARU") are assigned a point count method type because they are legacy data that were likely transcribed through continuous listening without the use of a spectrogram.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(wildrtrax) #to tidy data from wildtrax
library(readxl)

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#4. Set the WT version ----
v.wt <- "2026-06-02"

#5. Get the downloaded data object ----
load(file.path(root, "WildTrax", v.wt, paste0("01_wildtrax_raw_", v.wt, ".Rdata")))

bad_tasks = read_xlsx(file.path(root, "Dataset Assessment", "Exclusion", "Retenu_visite.xlsx"))

# HELPER FUNCTIONS ########

remove_bad_dates = function(dat, col = "date_time", begin_date = "1900-01-01", end_date = Sys.time()) {
  dat[dat[, col] >= begin_date & dat[, col] <= end_date, ]
}

remove_bad_coordinates = function(dat, xcol = "longitude", ycol = "latitude", lon_keep = c(-180, 180), lat_keep = c(-90, 90), lon_flip = c(180, 180), lat_flip = c(90, 90)) {
  xcol_values = dat[, xcol, drop = TRUE]
  ycol_values = dat[, ycol, drop = TRUE]
  
  dat = dat[!is.na(xcol_values) & !is.na(ycol_values), ]
  
  xcol_values_new = dat[, xcol, drop = TRUE]
  ycol_values_new = dat[, ycol, drop = TRUE]
  xcol_flip_bool = (xcol_values_new > lon_flip[1]) & (xcol_values_new < lon_flip[2])
  ycol_flip_bool = (ycol_values_new > lat_flip[1]) & (ycol_values_new < lat_flip[2])
  
  dat[xcol_flip_bool, xcol] = -dat[xcol_flip_bool, xcol]
  dat[ycol_flip_bool, ycol] = -dat[ycol_flip_bool, ycol]
  
  xcol_values_newer = dat[, xcol, drop = TRUE]
  ycol_values_newer = dat[, ycol, drop = TRUE]
  lat_lon_diff_values = abs(xcol_values_newer) != abs(ycol_values_newer)
  lat_lon_keep_values = (xcol_values_newer > lon_keep[1]) & (xcol_values_newer < lon_keep[2]) & (ycol_values_newer > lat_keep[1]) & (ycol_values_newer < lat_keep[2])
  
  dat[lat_lon_diff_values & lat_lon_keep_values, ]
}

#TIDY ARU DATA###########

#1. Collapse to single dataframe ----
aru <- do.call(rbind, aru.wt)

#2. Remove erroneous observations ----
# Locations: 
#           - some latitudes are positive when they should be negative -> swap those EXCEPT a few point counts from the Aleutians which should stay the same! AND there are some large negative longitudes (around the same area) that should be flipped! But not always reliable here so for now going to just leave out all large western longitudes. Suggest 171 W as a boundary and can look into it more later. Filter accordingly
#           - some longitudes are just the negative versions of the latitudes -> remove these (more generally remove when abs(lon) == abs(lat))
#           - remove anything with a location buffer (mostly ARU data) - these are incorrect locations
#           - some latitudes appear to be divided by 10 - don't want to fix that now but should be taken into account later

# Dates:
#           - remove any dates after the date in which the data were retrieved
#           - remove any dates in the year 1900 or earlier - oldest actual data appears to be BBS data from the 60's

# Other data qualities:
#           - remove anything where "task_is_complete" is not TRUE (only for ARU)
#           - remove any ARU detections after 10 minutes (and shorten survey duration to 10 minutes)
#           - remove duplicate instances of tag_id (for ARU) or survey_id-detection_distance-detection_time-species_code combinations (for PC)
#           - remove anything with more than 99.9% quantile of count for point count data for that species (or 10 if that quantile is smaller than 10)
MAX_ARU_TIME = 15 * 60
BEGIN_DATE = "1901-01-01"

aru.good = aru %>%
  # remove non-birds
  wt_tidy_species(remove=c("abiotic", "insect", "human")) %>%
  dplyr::filter(!(species_code %in% c("NONE")), !is.na(species_code)) %>%
  # estimate counts in the event of "too many to track" detections
  wt_replace_tmtt() %>%
  # remove any non-numeric values for individual_count
  dplyr::filter(individual_count > 0) %>%
  # remove erroneous noise
  dplyr::filter(max_noise_volume!="Extreme" | is.na(max_noise_volume), !max_noise_type %in% c("ARU Malfunction") | is.na(max_noise_type)) %>%
  # remove tasks labeled as bad by the "bad_tasks" dataframe
  left_join(bad_tasks, by = join_by(task_id == task_id, project_id == project_id, location == location, recording_date_time == recording_date_time)) %>%
  mutate(Retenu_Visite = ifelse(is.na(Retenu_Visite), "oui", Retenu_Visite)) %>%
  dplyr::filter(Retenu_Visite == "oui") %>%
  # remove anything with buffered locations or for which the task has not been completed yet
  dplyr::filter(is.na(location_buffer_m) | location_buffer_m == 0, task_is_complete %in% c("TRUE", "t")) %>%
  remove_bad_dates(col = "recording_date_time", begin_date = BEGIN_DATE, end_date = v.wt) %>%
  remove_bad_coordinates(lon_keep = c(-171, -52), lat_keep = c(30, 90), lon_flip = c(0, 170)) %>%
  # truncate detections to 15 minutes
  dplyr::filter(detection_time <= MAX_ARU_TIME) %>%
  mutate(task_duration = pmin(task_duration, MAX_ARU_TIME)) %>%
  # remove duplicate instances of the same tag_id
  distinct(tag_id, .keep_all = TRUE) %>%
  # remove duplicate detections of the same individual by grouping along "individual_order" and selecting the minimum (first) detection
  group_by(project_id, location_id, longitude, latitude, task_id, recording_id, recording_date_time, species_code, individual_order) %>%
  dplyr::filter(detection_time == min(detection_time)) %>%
  ungroup

#3. Tidy and format ----
#we have to filter to the first detection for each "individual_order" because some individuals have multiple tags
aru.tidy <- aru.good |> 
  rename(date_time = recording_date_time,
         duration = task_duration,
         method = task_method,
         survey_id = task_id,
         status = task_is_complete) |> 
  mutate(distance = Inf,
         sensor = "ARU",
         species = ifelse(species_code=="species", "UNKN", species_code)) |> 
  group_by(organization, project_id, location_id, location_buffer_m, longitude, latitude, survey_id, date_time, status, method, duration, distance, max_noise_type, max_noise_volume, species) |> 
  summarize(count = sum(individual_count)) |> 
  ungroup()

#TIDY PC DATA############

#1. Collapse to a single dataframe ----
pc <- do.call(rbind, pc.wt)

#2. Remove erroneous observations ----
MIN_PC_COUNT_CUTOFF = 10

pc.good = pc %>%
  # remove non-birds
  wt_tidy_species(remove=c("abiotic", "insect", "human")) %>%
  dplyr::filter(!(species_code %in% c("NONE")), !is.na(species_code)) %>%
  # remove tasks labeled as bad by the "bad_tasks" dataframe
  left_join(bad_tasks, by = join_by(survey_id == task_id, project_id == project_id, location == location, survey_date == recording_date_time)) %>%
  mutate(Retenu_Visite = ifelse(is.na(Retenu_Visite), "oui", Retenu_Visite)) %>%
  dplyr::filter(Retenu_Visite == "oui") %>%
  # remove anything with buffered locations or for which the task has not been completed yet
  dplyr::filter(is.na(location_buffer_m) | location_buffer_m == 0) %>%
  remove_bad_dates(col = "survey_date", begin_date = BEGIN_DATE, end_date = v.wt) %>%
  remove_bad_coordinates(lon_keep = c(-171, -52), lat_keep = c(30, 90), lon_flip = c(0, 170)) %>%
  # remove any surveys with an individual count of 0
  mutate(individual_count = as.numeric(individual_count),
         individual_count = ifelse(is.na(individual_count), 0, individual_count)) %>%
  group_by(survey_id) %>%
  dplyr::filter(all(individual_count > 0)) %>%
  ungroup %>%
  group_by(organization, project_id, location_id, longitude, latitude, survey_id, species_code) %>%
  mutate(total_count = sum(individual_count)) %>%
  ungroup

# remove any counts above the 99.9% quantile for the species (or 10)
pc.species.count.info = pc.good %>%
  group_by(organization, project_id, location_id, longitude, latitude, survey_id, survey_date, species_code) %>%
  mutate(total_count = sum(individual_count)) %>%
  ungroup %>% 
  distinct(organization, project_id, location_id, longitude, latitude, survey_id, survey_date, species_code, total_count) %>%
  group_by(species_code) %>%
  mutate(upper_q = pmax(MIN_PC_COUNT_CUTOFF, quantile(total_count, 0.999))) %>%
  ungroup %>%
  distinct(species_code, upper_q)

pc.good.final = pc.good %>%
  left_join(pc.species.count.info, by = join_by(species_code == species_code)) %>%
  dplyr::filter(total_count <= upper_q)

#3. Tidy and format ----
pc.tidy <- pc.good.final |> 
  wt_tidy_species(remove=c("abiotic", "insect", "human")) |> 
  rename(date_time = survey_date) |> 
  mutate(method = "PC",
         duration = as.integer(str_extract(survey_duration_method,
                                           "(?<=-)[0-9]+(?=min?)"))*60,
         distance = ifelse(str_sub(survey_distance_method, -3, -1) %in% c("INF", "ARU") | survey_distance_method=="UNKNOWN", Inf,
                           as.integer(str_extract(survey_distance_method,
                                                  "(?<=-)[0-9]+(?=m?)"))),
         max_noise_type = NA,
         max_noise_volume = NA,
         count = as.integer(individual_count),
         species = ifelse(species_code=="species", "UNKN", species_code),
         status = TRUE) |> 
  dplyr::select(all_of(colnames(aru.tidy)))

#PUT TOGETHER#########

#1. Combine ----
wt.tidy <- rbind(aru.tidy, pc.tidy)

#2. Filter out tasks we don't want ----

#filter to approximately North America
#clean up some bird codes
#only use species with 4 letter codes
#remove QC tasks that should not be used
#remove ARU tasks with too much noise
wt.use <- wt.tidy  |> 
  dplyr::filter(method != "None",
                !is.na(duration),
                !is.na(distance),
                !is.na(latitude),
                !is.na(date_time),
                str_length(species)==4,
                species!="4794") |> 
  mutate(species = case_when(species=="GRAJ" ~ "CAJA",
                             species=="PSFL" ~ "WEFL",
                             species=="MEGU" ~ "COGU",
                             !is.na(species) ~ species))
rm(wt.tidy)

# PRT : experimenting with some filters for removing errnoenous species reports, but going to leave this uncommented for now because takes a long time
# library(terra)
# wt_species_kdes = with(wt.use, make_individual_kdes(longitude, latitude, date_time, id = species, levels = 0.99, verbose = 1, min_locs = 50))
# 
# wt.use = wt.use %>%
#   arrange(species)
# wt_species_vect = vect(wt.use, geom = c("longitude", "latitude"), crs = add_EPSG(4326))
# 
# wt.use$d_kde_bound = unlist(lapply(unique(wt.use$species), function(sp) {
#   
#   message("beginning ", sp)
#   this_species = wt_species_vect[wt.use$species == sp, ]
#   this_kde = wt_species_kdes[wt_species_kdes$name == sp, ]
#   if (nrow(this_kde) == 0) return(rep(NA, nrow(this_species)))
#   
#   terra::distance(this_species, this_kde)
#   
# }))
# 
# d_nearest_conspecific = sapply(1:nrow(wt.use), function(row){
#   
#   if (row %% 100 == 0) message("beginning row ", row)
#   this_ind = wt_species_vect[row, ]
#   other_conspecifics = wt_species_vect[wt.use$species == this_ind$species, ]
#   conspecific_dists = terra::distance(other_conspecifics, this_ind)
#   
#   min(conspecific_dists[conspecific_dists > 0])
#   
# })

#3. Make wide ----
#we don't use wt_make_wide() because we're using a different format now
#remove columns we dont need anymore
wt.wide <- wt.use |> 
  pivot_wider(names_from=species, values_from=count, values_fn=sum, values_fill=0) |> 
  dplyr::select(-status, -location_buffer_m, -max_noise_type, -max_noise_volume)

#4. Save ----
save(wt.wide, aru.good, pc.good.final, file = file.path(root, "WildTrax", v.wt, paste0("02_wildtrax_clean_", v.wt, ".Rdata")))
