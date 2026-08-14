# ---
# title: BAM dataset - harmonize data
# author: Elly Knight
# created: March 4, 2026
# ---

#NOTES################################

#PURPOSE: This script wrangles & combines eBird data with the cleaned WildTrax data. Loading the ebird datasets takes time, be patient.
#TO DO: Make this flag instead of remove, and prioritize non-duplicates when flagging
#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(auk) #ebird data read
library(wildrtrax) #species list
library(data.table) #binding with missing columns

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Set the WildTrax version ----
v.wt <- "2026-07-10"

#4. Set the eBird version ----
v.ebd <- "Jun-2026"

#5. Get the WT data ----
load(file.path(
  root,
  "WildTrax",
  v.wt,
  paste0("02_wildtrax_clean_", v.wt, ".Rdata")
))

#6. Get the ebird data ----
ebd.ca <- read_ebd(file.path(
  root,
  "eBird",
  v.ebd,
  paste0("03_ebd_filtered_Canada_", v.ebd, ".txt")
))
ebd.mx <- read_ebd(file.path(
  root,
  "eBird",
  v.ebd,
  paste0("03_ebd_filtered_Mexico_", v.ebd, ".txt")
))
ebd.us <- read_ebd(file.path(
  root,
  "eBird",
  v.ebd,
  paste0("03_ebd_filtered_USA_", v.ebd, ".txt")
))
ebd.raw <- rbind(ebd.ca, ebd.mx, ebd.us)
rm(ebd.ca, ebd.mx, ebd.us)

#7. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#8. Set the column names ----
colnms <- colnames(wt.wide[, c(1:10)])

#9. Get the potential QPAD list ----
spp_qpad <- read.csv(file.path(root, "qpad_eligible_species_2026-07-10.csv"))

#HARMONIZE##############

#1. Species lookup ----
#take out duplicates of scientific name

dup <- c("GRAJ", "CORBRA", "MEGU", "PICHUD", "ANSROS", "PSFL")

spp_use <- wildrtrax::wt_get_species() |>
  dplyr::filter(
    species_class == "AVES",
    species_scientific_name != " ",
    !species_code %in% dup,
    species_code %in% spp_qpad$species_code
  ) |>
  rename(scientific_name = species_scientific_name) |>
  dplyr::select(species_code, scientific_name)

#2. Get unique checklists only ----
ebd.unique <- auk_unique(ebd.raw)
#rm(ebd.raw)

#5. Get surveys ----
ebd.visit <- ebd.unique |>
  dplyr::select(-c(common_name, scientific_name, observation_count)) |>
  unique()

#3. Tidy ebird data----
#Note this assumes observations with "X" individuals are 1s
#Filter out hotspots
#Replace common name with alpha code
#Filter to unique checklists only
#Filter to just the wt species so the object isn't so huge (and binds properly), make sure to zero fill
ebd.tidy <- ebd.unique |>
  dplyr::filter(
    locality_type != "H",
    !is.na(duration_minutes),
    !is.na(latitude)
  ) |>
  left_join(spp_wt) |>
  dplyr::filter(species_code %in% spp_use$species_code) |>
  right_join(ebd.visit) |>
  mutate(species_code = ifelse(is.na(species_code), "NONE", species_code)) |>
  mutate(
    source = "eBird",
    organization = "eBird",
    project_id = 99999,
    sensor = "PC",
    method = "eBird",
    buffer = 0,
    date_time = ymd_hms(
      paste0(observation_date, time_observations_started),
      tz = "America/Edmonton"
    ),
    distance = Inf,
    abundance = as.numeric(ifelse(
      observation_count == "X",
      1,
      observation_count
    )),
    duration = duration_minutes * 60
  ) |>
  rename(observer = observer_id, survey_id = checklist_id) |>
  dplyr::filter(!is.na(date_time))
rm(ebd.unique)

#4. Get location ids ----
ebd.loc <- ebd.visit |>
  dplyr::select(latitude, longitude) |>
  unique() |>
  mutate(location_id = row_number() + max(wt.wide$location_id))

#5. Make wide ----
ebd.wide <- ebd.tidy |>
  left_join(ebd.loc) |>
  dplyr::select(all_of(colnms), species_code, abundance) |>
  pivot_wider(
    names_from = species_code,
    values_from = abundance,
    values_fn = sum,
    values_fill = 0
  ) |>
  dplyr::select(-NONE)

#6. Put together ----
#sort the species
DEG_DECIMALS <- 3 # for rounding lat's and lon's - maybe if I want to be more precise I can eventually convert to UTMs and use a real spatial unit but I think this is fine for now
TIME_ROUND <- "1 minute"

all.wide <- rbindlist(list(wt.wide, ebd.wide), fill = TRUE)

all_wide_no_dups <- all.wide %>%
  mutate(
    lon_rounded = round(longitude, DEG_DECIMALS),
    lat_rounded = round(latitude, DEG_DECIMALS),
    time_rounded = round_date(date_time, TIME_ROUND),
    method_sort = as.integer(factor(
      method,
      levels = c("PC", "1SPT", "1SPM", "1SPM Audio/Visual hybrid", "eBird")
    ))
  ) %>%
  # sort by priority for keeping based on project
  # TO DO: Make this flag instead of remove, and prioritize non-duplicates when flagging
  arrange(method_sort, date_time) %>%
  group_by(lat_rounded, lon_rounded, time_rounded) %>%
  mutate(n_dups = n(), keep = c(1, numeric(n() - 1))) %>%
  ungroup %>%
  dplyr::filter(keep == 1)

dat <- all_wide_no_dups |>
  select(all_of(colnms), sort(setdiff(names(all.wide), all_of(colnms)))) |>
  mutate(across(-colnms, replace_na, 0))

#7. Save ----
save(
  dat,
  file = file.path(
    root,
    paste0("04_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, ".Rdata")
  )
)
