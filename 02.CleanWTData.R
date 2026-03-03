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

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData"

#3. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#4. Get the downloaded data object ----
v <- "2026-03-02"
load(file.path(root, "WildTrax", paste0("wildtrax_raw_", v, ".Rdata")))

#TIDY ARU DATA###########

#1. Collapse to single dataframe ----
aru <- do.call(rbind, aru.wt)

#2. Tidy and format ----
#we have to filter to the first detection for each "individual_order" because some individuals have multiple tags
aru.tidy <- aru |> 
  wt_tidy_species(remove=c("abiotic", "insect", "human")) |> 
  wt_replace_tmtt() |> 
  rename(date_time = recording_date_time,
         duration = task_duration,
         method = task_method,
         survey_id = task_id,
         status = task_is_complete) |> 
  mutate(distance = Inf,
         sensor = "ARU",
         species = ifelse(species_code=="species", "UNKN", species_code)) |> 
  group_by(organization, project_id, location_id, location_buffer_m, longitude, latitude, survey_id, date_time, status, method, duration, distance, max_noise_type, max_noise_volume, species, individual_order) |>
  dplyr::filter(detection_time==min(detection_time)) |> 
  group_by(organization, project_id, location_id, location_buffer_m, longitude, latitude, survey_id, date_time, status, method, duration, distance, max_noise_type, max_noise_volume, species) |> 
  summarize(count = sum(individual_count)) |> 
  ungroup()
  
#TIDY PC DATA############

#1. Collapse to a single dataframe ----
pc <- do.call(rbind, pc.wt)

#2. Tidy and format ----
pc.tidy <- pc |> 
  wt_tidy_species(remove=c("abiotic", "insect", "human")) |> 
  rename(date_time = survey_date) |> 
  mutate(method = "PC",
         duration = as.integer(str_extract(survey_duration_method,
                                           "(?<=-)[0-9]+(?=min?)"))*60,
         distance = ifelse(str_sub(survey_distance_method, -3, -1) %in% c("INF", "ARU"), Inf,
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
all.tidy <- rbind(aru.tidy, pc.tidy)

#2. Filter out tasks we don't want ----
#filter to approximately North America
all.use <- all.tidy  |> 
  dplyr::filter(method!="None",
                status %in% c("t", "TRUE"),
                (max_noise_volume!="Extreme" | is.na(max_noise_volume)),
                (!max_noise_type %in% c("ARU Malfunction") | is.na(max_noise_type)),
                (is.na(location_buffer_m) | location_buffer_m==0),
                !is.na(duration),
                !is.na(distance),
                !is.na(latitude),
                !is.na(date_time),
                latitude > 10,
                latitude < 85,
                longitude < -52,
                longitude > -168) 

#3. Make wide ----
#we don't use wt_make_wide() because we're using a different format now
#remove columns we dont need anymore
all.wide <- all.use |> 
  pivot_wider(names_from=species, values_from=count, values_fn=sum, values_fill=0) |> 
  dplyr::select(-status, -location_buffer_m, -max_noise_type, -max_noise_volume)

#4. Save ----
save(all.wide, file=file.path(root, "WildTrax", paste0("02_wildtrax_clean_", v, ".Rdata")))
