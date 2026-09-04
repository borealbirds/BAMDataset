# ---
# title: BAM dataset - harmonize data
# author: Elly Knight
# created: March 4, 2026
# ---

#NOTES################################

#PURPOSE: This script wrangles & combines eBird data with the cleaned WildTrax data. Loading the ebird datasets takes time, be patient.
#TO DO: Make this flag instead of remove, and prioritize non-duplicates when flagging

# PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(auk) #ebird data read
library(wildrtrax) #species list
library(data.table) #binding with missing columns
library(readxl)
library(terra)
library(sf)
library(duckdb)
library(duckspatial)

lapply(list.files("R", pattern = "\\.R", full.names = TRUE), source)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Set the WildTrax version ----
v.wt <- "2026-08-26"

#4. Set the eBird version ----
v.ebd <- "Jun-2026"

#5. Login to WildTrax----
source("WTlogin.R")
wt_auth()

# HELPER FUNCTIONS ########

remove_bad_dates <- function(
  dat,
  col = "date_time",
  begin_date = "1900-01-01",
  end_date = Sys.time()
) {
  dat[dat[, col] >= begin_date & dat[, col] <= end_date, ]
}

remove_bad_coordinates <- function(
  dat,
  xcol = "longitude",
  ycol = "latitude",
  lon_keep = c(-180, 180),
  lat_keep = c(-90, 90),
  lon_flip = c(180, 180),
  lat_flip = c(90, 90),
  flag = FALSE
) {
  if (flag && "flag_err_loc" %notin% names(dat)) {
    dat$flag_err_loc <- FALSE
  }

  xcol_values <- dat[, xcol, drop = TRUE]
  ycol_values <- dat[, ycol, drop = TRUE]

  dat <- dat[!is.na(xcol_values) & !is.na(ycol_values), ]

  xcol_values_new <- dat[, xcol, drop = TRUE]
  ycol_values_new <- dat[, ycol, drop = TRUE]
  xcol_flip_bool <- (xcol_values_new > lon_flip[1]) &
    (xcol_values_new < lon_flip[2])
  ycol_flip_bool <- (ycol_values_new > lat_flip[1]) &
    (ycol_values_new < lat_flip[2])

  dat[xcol_flip_bool, xcol] <- -dat[xcol_flip_bool, xcol]
  dat[ycol_flip_bool, ycol] <- -dat[ycol_flip_bool, ycol]

  xcol_values_newer <- dat[, xcol, drop = TRUE]
  ycol_values_newer <- dat[, ycol, drop = TRUE]
  lat_lon_diff_values <- abs(xcol_values_new) != abs(ycol_values_new)
  lat_lon_keep_values <- (xcol_values_new > lon_keep[1]) &
    (xcol_values_new < lon_keep[2]) &
    (ycol_values_new > lat_keep[1]) &
    (ycol_values_new < lat_keep[2])

  if (flag) {
    dat <- dat %>%
      mutate(
        flag_err_loc = flag_err_loc |
          !lat_lon_diff_values |
          !lat_lon_keep_values
      )
    return(dat)
  }

  dat[lat_lon_diff_values & lat_lon_keep_values, ]
}

wt_replace_tmtt <- function(data, calc = "round") {
  if (!"recording_date_time" %in% colnames(data)) {
    stop(
      "The `wt_replace_tmtt` function only works on data from the ARU sensor"
    )
  }
  check_none <- pull(distinct(select(data, species_code)))
  if (length(check_none) == 1 && check_none == "NONE") {
    stop("There are no species in this project")
  }
  .tmtt <- readRDS(system.file(
    "extdata",
    "tmtt_predictions.rds",
    package = "wildrtrax"
  ))
  dat.tmtt <- mutate(data, id = row_number())
  dat.tmt <- filter(dat.tmtt, abundance %in% c("TMTT", "TNPE"))
  if (nrow(dat.tmt) > 0) {
    dat.tmt <- select(
      mutate(
        inner_join(
          mutate(
            dat.tmt,
            species_code = ifelse(
              species_code %in% .tmtt$species_code,
              species_code,
              "species"
            ),
            observer_id = as.integer(ifelse(
              observer_id %in%
                .tmtt$observer_id,
              observer_id,
              0
            ))
          ),
          select(.tmtt, species_code, observer_id, pred),
          by = c("species_code", "observer_id")
        ),
        abundance = case_when(
          calc == "round" ~ round(pred),
          calc == "ceiling" ~ ceiling(pred),
          calc == "floor" ~ floor(pred),
          TRUE ~ NA_real_
        )
      ),
      -pred
    )
  }

  dat.tmtt <- dat.tmtt %>%
    mutate(
      abundance = case_when(
        abundance %in% c("TMTT", "TNPE") ~ NA_real_,
        .default = as.numeric(abundance)
      )
    ) %>%
    dplyr::select(-id)

  return(dat.tmtt)
}

# WILDTRAX DATA WRANGLING #################

#1. Get the downloaded data object ----
db_conn_path <- file.path(
  root,
  "WildTrax",
  v.wt,
  paste0("01_wildtrax_raw_", v.wt, ".duckdb")
)
db_conn <- dbConnect(duckdb(db_conn_path))

bad_tasks <- read_xlsx(file.path(
  root,
  "Dataset Assessment",
  "Exclusion",
  "Retenu_visite.xlsx"
))
bad_aru_projects <- read_csv(file.path(
  root,
  "Dataset Assessment",
  "Exclusion",
  "aru_bad_timestamps_prt_20260716.csv"
))

#2. Tidy ARU data ----
#we have to filter to the first detection for each "individual_order" because some individuals have multiple tags
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
#           - flag counts above a sample-size-appropriate species quantile, with a minimum count cutoff of 10
MAX_ARU_TIME <- 10 * 60
BEGIN_DATE <- "1901-01-01"
PERC_TOD_THRESHOLD <- 0.01 # we remove any ARU counts that come from a time of day that represents less than this percentage of possible time of day bins
MIN_COUNT_CUTOFF <- 10
MIN_Q995_OBSERVATIONS <- 200
MIN_Q999_OBSERVATIONS <- 1000

timeofday_cov <- cov_tod_bin("timeofday")

aru_qpad_ready <- ddbs_read_table(db_conn, "aru") %>%
  # add lat and lon as columns so they're not just in geometries
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude = st_coordinates(.)[, 2],
    location = as.character(location),
    recording_date_time = as.POSIXct(
      recording_date_time,
      format = "%Y-%m-%d %H:%M:%OS",
      tz = "UTC"
    )
  ) %>%
  # drop geometry for now
  st_drop_geometry %>%
  # remove non-birds
  wt_tidy_species(remove = c("abiotic", "insect", "human")) %>%
  dplyr::filter(
    !(species_code %in% c("NONE")),
    !is.na(species_code),
    !is.na(task_duration),
    # Records explicitly labelled with no survey method are not valid ARU surveys.
    is.na(task_method) |
      str_to_upper(str_squish(task_method)) != "NONE"
  ) %>%
  # estimate counts in the event of "too many to track" detections
  wt_replace_tmtt() |>
  # remove any non-numeric values for abundance
  dplyr::filter(abundance > 0) |>
  mutate(
    # flag sightings that have missing timestamp information (often shows up as being recorded at midnight)
    flag_err_date = hour(recording_date_time) == 0 &
      minute(recording_date_time) == 0,
    # flag sightings with erroneous noise
    flag_noise = (!is.na(max_noise_volume) & max_noise_volume == "Extreme") |
      (!is.na(max_noise_type) & max_noise_type == "ARU Malfunction")
  ) %>%
  # remove tasks labeled as bad by the "bad_tasks" dataframes
  left_join(
    bad_tasks,
    by = join_by(
      task_id == task_id,
      project_id == project_id,
      location == location,
      recording_date_time == recording_date_time
    )
  ) %>%
  mutate(
    Retenu_Visite = ifelse(is.na(Retenu_Visite), "oui", Retenu_Visite),
    # flag any tasks that have been manually recommended for removal
    flag_manual = Retenu_Visite != "oui" |
      project_id %in% bad_aru_projects$project_id,
    # flag anything with buffered locations
    flag_err_loc = !is.na(location_buffer_m) & location_buffer_m > 0,
    # flag tasks labeled as incomplete
    flag_task_wt = task_is_complete %notin% c("TRUE", "t")
  ) %>%
  # more flagging for locations - if coordinates make no sense
  remove_bad_coordinates(
    lon_keep = c(-171, -52),
    lat_keep = c(30, 90),
    lon_flip = c(0, 170),
    flag = TRUE
  ) %>%
  # truncate detections and survey duration to 10 minutes
  dplyr::filter(detection_time <= MAX_ARU_TIME) %>%
  mutate(task_duration = pmin(task_duration, MAX_ARU_TIME)) %>%
  # remove duplicate instances of the same tag_id
  distinct(tag_id, .keep_all = TRUE) %>%
  # remove duplicate detections of the same individual by grouping along "individual_order" and selecting the minimum (first) detection
  group_by(
    project_id,
    location_id,
    longitude,
    latitude,
    task_id,
    recording_id,
    recording_date_time,
    species_code,
    individual_order
  ) %>%
  dplyr::filter(detection_time == min(detection_time)) %>%
  ungroup %>%
  # rename task ID's so they don't overlap with point counts
  mutate(task_id = paste0("ARU_", task_id)) %>%
  # remove "outlier" counts that take place at the wrong time of day. this involves calculating time of day which is a bit tedious but doable
  mutate(
    timeofday = timeofday_cov$get(longitude, latitude, recording_date_time)
  ) %>%
  group_by(project_id, timeofday) %>%
  mutate(n_tod = n()) %>%
  ungroup %>%
  group_by(project_id) %>%
  mutate(n_proj = n()) %>%
  ungroup %>%
  mutate(
    p_tod = n_tod / n_proj,
    flag_err_date = flag_err_date | p_tod < PERC_TOD_THRESHOLD
  ) |>
  rename(
    date_time = recording_date_time,
    duration = task_duration,
    method = task_method,
    survey_id = task_id,
    status = task_is_complete
  ) |>
  mutate(
    data_source = "WildTrax",
    sensor = "ARU",
    distance = Inf,
    species = ifelse(species_code == "species", "UNKN", species_code)
  ) %>%
  collect

aru.tidy <- aru_qpad_ready |>
  group_by(
    data_source,
    sensor,
    organization,
    project_id,
    location_id,
    location_buffer_m,
    longitude,
    latitude,
    survey_id,
    date_time,
    status,
    method,
    duration,
    distance,
    max_noise_type,
    max_noise_volume,
    species,
    flag_err_date,
    flag_err_loc,
    flag_noise,
    flag_task_wt,
    flag_manual
  ) |>
  summarize(count = sum(abundance)) |>
  ungroup()

#3. Tidy point count data ----

pc_qpad_ready <- ddbs_read_table(db_conn, "pc") %>%
  # add lat and lon as columns so they're not just in geometries
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude = st_coordinates(.)[, 2],
    location = as.character(location),
    survey_date_time = as.POSIXct(
      survey_date_time,
      format = "%Y-%m-%d %H:%M:%OS",
      tz = "UTC"
    ),
    survey_date = survey_date_time
  ) %>%
  # drop geometry for now
  st_drop_geometry %>%
  # remove non-birds
  wt_tidy_species(remove = c("abiotic", "insect", "human")) %>%
  # remove surveys with no species or duration info (no distance info is OK, we just assume all birds were counted regardless of distance)
  dplyr::filter(
    !(species_code %in% c("NONE")),
    !is.na(species_code),
    survey_duration_method != "UNKNOWN"
  ) |>
  mutate(
    # flag sightings that have missing timestamp information (often shows up as being recorded at midnight)
    flag_err_date = hour(survey_date_time) == 0 & minute(survey_date_time) < 2
  ) %>%
  # remove tasks labeled as bad by the "bad_tasks" dataframe
  left_join(
    bad_tasks,
    by = join_by(
      survey_id == task_id,
      project_id == project_id,
      location == location,
      survey_date_time == recording_date_time
    )
  ) %>%
  mutate(
    Retenu_Visite = ifelse(is.na(Retenu_Visite), "oui", Retenu_Visite),
    # flag any tasks that have been manually recommended for removal
    flag_manual = Retenu_Visite != "oui",
    # flag anything with buffered locations
    flag_err_loc = !is.na(location_buffer_m) & location_buffer_m > 0
  ) %>%
  # more flagging for locations - if coordinates make no sense
  remove_bad_coordinates(
    lon_keep = c(-171, -52),
    lat_keep = c(30, 90),
    lon_flip = c(0, 170),
    flag = TRUE
  ) %>%
  # remove any surveys with an individual count of 0
  mutate(
    abundance = as.numeric(abundance),
    abundance = ifelse(is.na(abundance), 0, abundance)
  ) %>%
  group_by(survey_id) %>%
  dplyr::filter(all(abundance > 0)) %>%
  ungroup %>%
  # rename task ID's so they don't overlap with ARU surveys
  mutate(survey_id = paste0("PC_", survey_id)) |>
  rename(date_time = survey_date_time, abundance = abundance) %>%
  collect

pc.tidy <- pc_qpad_ready %>%
  group_by(
    organization,
    project_id,
    location_id,
    longitude,
    latitude,
    survey_id,
    species_code
  ) %>%
  mutate(total_count = sum(abundance)) %>%
  ungroup %>%
  mutate(
    data_source = "WildTrax",
    sensor = "PC",
    method = "PC",
    duration = as.integer(str_extract(
      survey_duration_method,
      "(?<=-)[0-9]+(?=min?)"
    )) *
      60,
    distance = ifelse(
      str_sub(survey_distance_method, -3, -1) %in%
        c("INF", "ARU") |
        survey_distance_method == "UNKNOWN",
      Inf,
      as.integer(str_extract(survey_distance_method, "(?<=-)[0-9]+(?=m?)"))
    ),
    max_noise_type = NA_character_,
    max_noise_volume = NA_character_,
    count = as.integer(abundance),
    species = ifelse(species_code == "species", "UNKN", species_code),
    status = TRUE
  ) %>%
  mutate(
    # these two are only for ARUs
    flag_noise = FALSE,
    flag_task_wt = FALSE
  ) |>
  dplyr::select(all_of(colnames(aru.tidy)))

#4. Combine ----
wt.tidy <- list(aru.tidy, pc.tidy) |>
  map(
    ~ .x |>
      mutate(
        across(
          all_of(c(
            "data_source",
            "sensor",
            "organization",
            "survey_id",
            "method",
            "max_noise_type",
            "max_noise_volume",
            "species"
          )),
          as.character
        ),
        across(
          all_of(c(
            "project_id",
            "location_id",
            "location_buffer_m",
            "longitude",
            "latitude",
            "duration",
            "distance",
            "count"
          )),
          as.double
        ),
        across(
          all_of(c(
            "status",
            "flag_err_date",
            "flag_err_loc",
            "flag_noise",
            "flag_task_wt",
            "flag_manual"
          )),
          as.logical
        )
      )
  ) |>
  bind_rows() |>
  # Retain valid four-character species codes and harmonize legacy WildTrax
  # codes before survey-species counts are combined.
  mutate(
    species = case_when(
      species == "GRAJ" ~ "CAJA",
      species == "PSFL" ~ "WEFL",
      species == "MEGU" ~ "COGU",
      TRUE ~ species
    )
  ) |>
  dplyr::filter(
    str_length(species) == 4,
    species != "4794"
  ) %>%
  # group together observations of the same species that are split up (e.g., from point counts with different distances) and sum
  group_by(survey_id, species) %>%
  mutate(count = sum(count)) %>%
  ungroup %>%
  distinct(survey_id, species, .keep_all = TRUE)

#EBIRD DATA WRANGLING #################

ebd.unique <- fread(file.path(
  root,
  "eBird",
  v.ebd,
  paste0("03_ebd_filtered_ALL_", v.ebd, ".csv")
))

spp_qpad <- read.csv(file.path(root, "qpad_eligible_species_2026-08-26.csv"))

#take out duplicates of scientific name
dup <- c("GRAJ", "CORBRA", "MEGU", "PICHUD", "ANSROS", "PSFL")

spp_use <- wildrtrax::wt_get_species() |>
  dplyr::filter(
    species_class == "AVES",
    !is.na(species_scientific_name),
    str_squish(species_scientific_name) != "",
    !species_code %in% dup,
    species_code %in% spp_qpad$species_code
  ) |>
  transmute(
    species_code,
    scientific_name_key = str_to_upper(str_squish(species_scientific_name))
  ) |>
  distinct()


#1. Tidy eBird data ----
# eBird times are local clock times, but the filtered files do not retain a
# timezone. UTC is used here as a neutral storage timezone; this does not
# convert the observations from local time to actual UTC.
# Observations reported as "X" are conservatively assigned a count of one.
ebd.tidy <- ebd.unique %>%
  mutate(
    number_observers = ifelse(
      "number_observers" %notin% names(.),
      1,
      number_observers
    )
  ) %>%
  dplyr::filter(
    # only include non-hotspots because hotspots are often defined at spatially imprecise locations
    locality_type != "H",
    # make sure nothing is NA or missing
    !is.na(checklist_id),
    checklist_id != "",
    !is.na(duration_minutes),
    duration_minutes > 0,
    !is.na(latitude),
    !is.na(longitude),
    # don't want checklists from large parties as these probably can't be compared to point counts
    number_observers == 1
  ) |>
  mutate(
    scientific_name_key = str_to_upper(str_squish(scientific_name)),
    date_time = ymd_hms(
      paste(observation_date, time_observations_started),
      tz = "UTC",
      quiet = TRUE
    ),
    count = case_when(
      observation_count == "X" ~ 1,
      TRUE ~ suppressWarnings(as.numeric(observation_count))
    )
  ) |>
  dplyr::filter(!is.na(date_time), !is.na(count), count > 0) |>
  left_join(
    spp_use,
    by = "scientific_name_key",
    relationship = "many-to-one"
  ) |>
  dplyr::filter(!is.na(species_code)) |>
  transmute(
    data_source = "eBird",
    sensor = "PC",
    organization = "eBird",
    project_id = 99999,
    location_id = suppressWarnings(
      max(c(0, wt.tidy$location_id), na.rm = TRUE)
    ) +
      dense_rank(paste(latitude, longitude, sep = ":")),
    location_buffer_m = 0,
    longitude,
    latitude,
    survey_id = paste0("EBIRD_", checklist_id),
    date_time,
    status = TRUE,
    method = "eBird",
    duration = as.numeric(duration_minutes) * 60,
    distance = Inf,
    max_noise_type = NA_character_,
    max_noise_volume = NA_character_,
    species = species_code,
    count,
    flag_err_date = FALSE,
    flag_err_loc = FALSE,
    flag_noise = FALSE,
    flag_task_wt = FALSE,
    flag_manual = FALSE
  ) |>
  dplyr::select(all_of(names(wt.tidy)))

#rm(ebd.unique, ebd.raw)

#COMBINE##############

#1. Put the three sources together and apply common flags ----
# The source versions provide reproducible upper date limits. This avoids
# changing flags merely because the script is rerun on a later date.
EBD_END_DATE <- as.Date(ceiling_date(my(v.ebd), "month") - days(1))

all.tidy <- bind_rows(wt.tidy, ebd.tidy) |>
  mutate(
    .source_end_date = case_when(
      data_source == "WildTrax" ~ as.Date(v.wt),
      data_source == "eBird" ~ EBD_END_DATE,
      TRUE ~ as.Date(NA)
    ),
    # Shared breeding-season flag: June 1 through July 15, inclusive.
    flag_doy = is.na(date_time) |
      yday(date_time) %notin% seq(152, 196),
    # Retain source-specific date flags and add common date validity checks.
    flag_err_date = replace_na(flag_err_date, FALSE) |
      is.na(date_time) |
      as.Date(date_time) < as.Date(BEGIN_DATE) |
      as.Date(date_time) > .source_end_date,
    # Retain source-specific location flags and add common physical checks.
    flag_err_loc = replace_na(flag_err_loc, FALSE) |
      is.na(longitude) |
      is.na(latitude) |
      longitude < -180 |
      longitude > 180 |
      latitude < -90 |
      latitude > 90 |
      abs(longitude) == abs(latitude),
    flag_noise = replace_na(flag_noise, FALSE),
    flag_task_wt = replace_na(flag_task_wt, FALSE),
    flag_manual = replace_na(flag_manual, FALSE)
  ) |>
  dplyr::select(-.source_end_date) |>
  mutate(
    # Preliminary survey-level flag. Duplicate status is added when visit is
    # created; species-level count flags are stored separately in bird.
    flag_any = if_any(
      all_of(c(
        "flag_doy",
        "flag_err_date",
        "flag_err_loc",
        "flag_noise",
        "flag_task_wt",
        "flag_manual"
      )),
      identity
    )
  )

# A survey must have exactly one data source and sensor assignment. Output object should have zero rows.
survey_provenance_issues <- all.tidy |>
  distinct(survey_id, data_source, sensor) |>
  count(survey_id, name = "n_provenance_combinations") |>
  dplyr::filter(n_provenance_combinations > 1)
nrow(survey_provenance_issues)

# The long-data key is one row per survey and species. Output object should have zero rows.
survey_species_issues <- all.tidy |>
  count(survey_id, species, name = "n_rows") |>
  dplyr::filter(n_rows > 1)
nrow(survey_species_issues)

#2. Split into visit and bird tables ----

# One row per survey, containing only survey-level metadata and flags. When an eBird checklist group has multiple metadata variants, retain the variant represented by the most species. The remaining fields provide deterministic tie-breakers.
visit_core <- all.tidy |>
  dplyr::select(-count) |>
  group_by(across(-species)) |>
  summarize(n_species = n(), .groups = "drop") |>
  arrange(
    survey_id,
    desc(n_species),
    date_time,
    latitude,
    longitude,
    duration,
    location_id
  ) |>
  distinct(survey_id, .keep_all = TRUE) |>
  dplyr::select(-n_species) |>
  # Candidate duplicate visits have approximately the same location and start time. Only groups containing WildTrax are flagged, which avoids treating independent eBird-only checklists as duplicates.
  mutate(
    .dup_latitude = round(latitude, 3),
    .dup_longitude = round(longitude, 3),
    .dup_date_time = round_date(date_time, "10 minutes")
  ) |>
  # Retain one survey per duplicate group by prioritizing: (1) surveys with no existing flags, (2) WildTrax over eBird, (3) ARU over point count, (4) longer duration, (5) greater survey distance, and (6) survey_id as a deterministic final tie-breaker. All other surveys in the group are flagged.
  arrange(
    .dup_date_time,
    .dup_latitude,
    .dup_longitude,
    flag_any,
    desc(data_source == "WildTrax"),
    desc(sensor == "ARU"),
    desc(duration),
    desc(distance),
    survey_id
  ) |>
  group_by(.dup_latitude, .dup_longitude, .dup_date_time) |>
  mutate(
    flag_dup = n() > 1 &
      any(data_source == "WildTrax") &
      row_number() > 1,
    # flag_any remains a survey-level flag: filtering it removes the entire
    # visit, including secondary copies from duplicate groups.
    flag_any = flag_any | flag_dup
  ) |>
  ungroup() |>
  arrange(survey_id) |>
  dplyr::select(-starts_with(".dup_"))

# Mark visits within common study regions. Perform each spatial predicate once
# per distinct location, rather than once per visit, and transform the points
# to each region's native CRS. The model-region layer contains separate BCR
# polygons, so dissolve it before testing point membership.
model_region <- st_read(
  file.path(root, "regions", "model", "BAM_BCRNMv5_3978.shp"),
  quiet = TRUE
) |>
  st_geometry() |>
  st_union()

canada_region <- st_read(
  file.path(root, "regions", "Canada", "CAN_adm0.shp"),
  quiet = TRUE
) |>
  st_geometry()

repair_invalid <- function(region) {
  geometry_is_valid <- st_is_valid(region)
  rbind(
    region[geometry_is_valid, ],
    st_make_valid(region[!geometry_is_valid, ])
  )
}

boreal_region <- st_read(
  file.path(root, "regions", "boreal", "NABoreal.shp"),
  quiet = TRUE
) |>
  repair_invalid() |>
  st_geometry() |>
  st_union()

visit_locations_sf <- visit_core |>
  dplyr::distinct(longitude, latitude) |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

inside_region <- function(points, region) {
  transformed_points <- st_transform(points, st_crs(region))
  lengths(st_intersects(transformed_points, region, prepared = TRUE)) > 0
}

visit_locations <- visit_locations_sf |>
  mutate(
    in_modelregion = inside_region(visit_locations_sf, model_region),
    in_Canada = inside_region(visit_locations_sf, canada_region),
    in_boreal = inside_region(visit_locations_sf, boreal_region)
  ) |>
  st_drop_geometry()

visit <- visit_core |>
  left_join(
    visit_locations,
    by = c("longitude", "latitude")
  )

# One row per survey and species, containing only bird detections. The
# flag_err_count field applies only to that species detection and is therefore
# not included in visit$flag_any. Use the
# 99.9% quantile for species with at least 1,000 observations, the 99.5%
# quantile for species with 200-999 observations, and no automated count flag
# below 200 observations. Retain only the logical flag in the full table.
bird <- all.tidy |>
  dplyr::select(survey_id, species, count) |>
  group_by(species) |>
  mutate(
    flag_err_count = if (n() >= MIN_Q999_OBSERVATIONS) {
      count >
        max(
          MIN_COUNT_CUTOFF,
          quantile(count, probs = 0.999, na.rm = TRUE, names = FALSE)
        )
    } else if (n() >= MIN_Q995_OBSERVATIONS) {
      count >
        max(
          MIN_COUNT_CUTOFF,
          quantile(count, probs = 0.995, na.rm = TRUE, names = FALSE)
        )
    } else {
      rep(FALSE, n())
    }
  ) |>
  ungroup()

#3. Check & interrogate output ----

# Each survey should have exactly one row in visit. Output object should have zero rows.
visit_key_issues <- visit |>
  count(survey_id, name = "n_rows") |>
  dplyr::filter(n_rows > 1)
nrow(visit_key_issues)

# Every bird record should link to a survey in visit. Output object should have zero rows.
bird_visit_issues <- bird |>
  anti_join(
    dplyr::select(visit, survey_id),
    by = "survey_id"
  ) |>
  distinct(survey_id)
nrow(bird_visit_issues)

# Each bird detection should have a unique survey and species key. Output object should have zero rows.
bird_key_issues <- bird |>
  count(survey_id, species, name = "n_rows") |>
  dplyr::filter(n_rows > 1)
nrow(bird_key_issues)

# Visit flags should be complete and represented correctly by flag_any. Output
# object should have zero rows.
visit_flag_issues <- visit |>
  dplyr::filter(
    if_any(starts_with("flag_"), is.na) |
      flag_any !=
        if_any(
          all_of(c(
            "flag_doy",
            "flag_err_date",
            "flag_err_loc",
            "flag_noise",
            "flag_task_wt",
            "flag_manual",
            "flag_dup"
          )),
          identity
        )
  )
nrow(visit_flag_issues)

# Every visit should receive a value for each spatial membership field. Output
# object should have zero rows.
spatial_membership_issues <- visit |>
  dplyr::filter(
    if_any(
      all_of(c("in_modelregion", "in_Canada", "in_boreal")),
      is.na
    )
  )
nrow(spatial_membership_issues)

# Bird values and count flags should be complete and valid. Output object should have zero rows.
bird_value_issues <- bird |>
  dplyr::filter(
    is.na(species) |
      species == "" |
      is.na(count) |
      !is.finite(count) |
      count <= 0 |
      is.na(flag_err_count)
  )
nrow(bird_value_issues)

# Provide one compact overview; all checks should report zero issues.
output_check_summary <- tibble(
  check = c(
    "survey provenance",
    "pre-split survey-species key",
    "visit key",
    "bird key",
    "bird-to-visit link",
    "visit flags",
    "spatial membership",
    "bird values"
  ),
  n_issues = c(
    nrow(survey_provenance_issues),
    nrow(survey_species_issues),
    nrow(visit_key_issues),
    nrow(bird_key_issues),
    nrow(bird_visit_issues),
    nrow(visit_flag_issues),
    nrow(spatial_membership_issues),
    nrow(bird_value_issues)
  )
)
output_check_summary

# Retain a compact review table containing both the selected survey and the
# flagged copies from each duplicate group.
duplicate_visit_review <- visit |>
  mutate(
    .dup_latitude = round(latitude, 3),
    .dup_longitude = round(longitude, 3),
    .dup_date_time = round_date(date_time, "10 minutes")
  ) |>
  group_by(.dup_latitude, .dup_longitude, .dup_date_time) |>
  dplyr::filter(any(flag_dup)) |>
  mutate(duplicate_group = cur_group_id()) |>
  ungroup() |>
  dplyr::select(
    duplicate_group,
    survey_id,
    flag_dup,
    data_source,
    sensor,
    organization,
    project_id,
    location_id,
    date_time,
    latitude,
    longitude,
    duration,
    distance,
    flag_doy,
    flag_err_date,
    flag_err_loc,
    flag_noise,
    flag_task_wt,
    flag_manual,
    flag_any
  ) |>
  arrange(duplicate_group, flag_dup, survey_id)

# Show the number of retained and flagged surveys from each source and sensor.
duplicate_visit_summary <- visit |>
  count(data_source, sensor, flag_dup, name = "n_surveys")
duplicate_visit_summary

#4. Add a lookup table for truncation values ----

# Retain one compact diagnostic row per species rather than attaching these thresholds to every row in the full bird table.
count_limits <- bird |>
  group_by(species) |>
  summarize(
    n_observations = n(),
    quantile_probability = if (n_observations >= MIN_Q999_OBSERVATIONS) {
      0.999
    } else if (n_observations >= MIN_Q995_OBSERVATIONS) {
      0.995
    } else {
      NA_real_
    },
    count_quantile = if (n_observations >= MIN_Q999_OBSERVATIONS) {
      quantile(
        count,
        probs = 0.999,
        na.rm = TRUE,
        names = FALSE
      )
    } else if (n_observations >= MIN_Q995_OBSERVATIONS) {
      quantile(
        count,
        probs = 0.995,
        na.rm = TRUE,
        names = FALSE
      )
    } else {
      NA_real_
    },
    applied_cutoff = max(MIN_COUNT_CUTOFF, count_quantile),
    max_count = max(count, na.rm = TRUE),
    max_to_cutoff_ratio = max_count / applied_cutoff,
    n_flagged = sum(flag_err_count, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_flagged), desc(max_count))

#5. Save ----

UTM_CRS <- 4326

aru_sf <- st_as_sf(
  aru_qpad_ready,
  coords = c("longitude", "latitude"),
  crs = UTM_CRS
) %>%
  mutate(longitude = st_coordinates(.)[, 1], latitude = st_coordinates(.)[, 2])
pc_sf <- st_as_sf(
  pc_qpad_ready,
  coords = c("longitude", "latitude"),
  crs = UTM_CRS
) %>%
  mutate(longitude = st_coordinates(.)[, 1], latitude = st_coordinates(.)[, 2])
visit_sf <- st_as_sf(
  visit,
  coords = c("longitude", "latitude"),
  crs = UTM_CRS
) %>%
  mutate(longitude = st_coordinates(.)[, 1], latitude = st_coordinates(.)[, 2])

## write to duckdb databases - for QPAD and for models ----

db_conn_path_qpad <- file.path(
  root,
  paste0("03_QPADDataset_WT-", v.wt, "_EBd-", v.ebd, ".duckdb")
)
db_conn_path_models <- file.path(
  root,
  paste0("03_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, ".duckdb")
)

### qpad first! ----

# make connection
ddbs_write_dataset(aru_sf, db_conn_path_qpad, layer = "aru", overwrite = TRUE)

# re-load connection to be able to add more data
db_conn_qpad <- dbConnect(duckdb(db_conn_path_qpad))
dbExecute(db_conn_qpad, "LOAD spatial") # have to do this for whatever reason

# write the other file
ddbs_write_table(db_conn_qpad, pc_sf, name = "pc", overwrite = TRUE)

### now for the models ----

# make connection
ddbs_write_dataset(
  visit_sf,
  db_conn_path_models,
  layer = "visit",
  overwrite = TRUE
)

# re-load connection to be able to add more data
db_conn_models <- dbConnect(duckdb(db_conn_path_models))
dbExecute(db_conn_models, "LOAD spatial") # have to do this for whatever reason

# write the other files
dbWriteTable(db_conn_models, name = "bird", bird)
dbWriteTable(db_conn_models, name = "count_limits", count_limits)
dbWriteTable(db_conn_models, name = "output_check_summary", output_check_summary)
