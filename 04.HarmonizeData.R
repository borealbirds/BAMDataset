# ---
# title: BAM dataset - harmonize data
# author: Elly Knight
# created: March 4, 2026
# ---

#NOTES################################

#PURPOSE: This script wrangles combines eBird data with the cleaned WildTrax data.

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
v.ebd <- "Jan-2026"

#5. Get the data ----
load(file.path(root, "WildTrax", v.wt, paste0("02_wildtrax_clean_", v.wt, ".Rdata")))
ebd.raw <- read_ebd(file.path(root, "eBird", v.ebd, paste0("03_ebd_filtered_CA_", v.ebd, ".txt")))

#6. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#7. Set the column names ----
colnms <- colnames(wt.wide[,c(1:10)])

#HARMONIZE##############

#1. Species lookup ----
#take out duplicates of scientific name

dup <- c("GRAJ", "CORBRA", "MEGU", "PICHUD", "ANSROS", "PSFL")

spp_wt <- wildrtrax::wt_get_species() |> 
  dplyr::filter(species_class=="AVES",
                species_scientific_name!=" ",
                !species_code %in% dup) |> 
  rename(scientific_name = species_scientific_name) |> 
  dplyr::select(species_code, scientific_name)

#2. Get unique checklists only ----
ebd.unique <- auk_unique(ebd.raw)

#3. Tidy ebird data----
#Note this assumes observations with "X" individuals are 1s
#Filter out hotspots
#Replace common name with alpha code
#Filter to unique checklists only
#Species with 4 letter alpha codes only
ebd.tidy <- ebd.unique |> 
  mutate(source = "eBird",
         organization = "eBird",
         project_id=99999,
         sensor="PC",
         method="eBird",
         buffer=0,
         date_time = ymd_hms(paste0(observation_date, time_observations_started), tz="America/Edmonton"),
         distance = Inf,
         abundance = as.numeric(ifelse(observation_count=="X", 1, observation_count)),
         duration = duration_minutes*60) |> 
  rename(observer = observer_id,
         survey_id = checklist_id) |> 
  left_join(spp_wt) |> 
  dplyr::filter(locality_type!="H",
                !is.na(duration),
                !is.na(distance),
                !is.na(latitude),
                !is.na(date_time),
                str_length(species_code)==4)

#4. Get location ids ----  
ebd.loc <- ebd.tidy |> 
  dplyr::select(latitude, longitude) |> 
  unique() |> 
  mutate(location_id = row_number() + max(wt.wide$location_id))
  
#5. Make wide ----
ebd.wide <- ebd.tidy |> 
  left_join(ebd.loc) |>
  dplyr::select(all_of(colnms), species_code, abundance) |> 
  pivot_wider(names_from=species_code, values_from=abundance, values_fn=sum, values_fill=0)

#6. Put together ----
#sort the species
DEG_DECIMALS = 3 # for rounding lat's and lon's - maybe if I want to be more precise I can eventually convert to UTMs and use a real spatial unit but I think this is fine for now
TIME_ROUND = "1 minute"

all.wide <- rbindlist(list(wt.wide, ebd.wide), fill=TRUE)

all_wide_no_dups = all.wide %>%
  mutate(lon_rounded = round(longitude, DEG_DECIMALS),
         lat_rounded = round(latitude, DEG_DECIMALS),
         time_rounded = round_date(date_time, TIME_ROUND),
         method_sort = as.integer(factor(method, levels = c("PC", "1SPT", "1SPM", "1SPM Audio/Visual hybrid", "eBird")))) %>%
  # sort by priority for keeping based on project
  arrange(method_sort, date_time) %>%
  group_by(lat_rounded, lon_rounded, time_rounded) %>%
  mutate(n_dups = n(),
         keep = c(1, numeric(n() - 1))) %>%
  ungroup %>%
  dplyr::filter(keep == 1)

dat <- all_wide_no_dups |> 
  select(all_of(colnms), sort(setdiff(names(all.wide), all_of(colnms)))) |> 
  mutate(across(-colnms, replace_na, 0))

#7. Save ----
save(dat, file=file.path(root, paste0("04_BAMDataset_WT-", v.wt, "_EBd-", v.ebd,  ".Rdata")))
