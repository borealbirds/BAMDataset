# ---
# title: BAM dataset - filter eBird data
# author: Elly Knight
# created: March 4, 2026
# ---

#NOTES################################

#PURPOSE: This script wrangles raw eBird datasets for combination with wildtrax data.

#raw eBird data is downloaded from the eBird interface at https://ebird.org/data/download/ebd prior to wrangling with the "auk" package and will require a request for access. Use the custom download tool to download only the datasets for Canada and the US instead of the global dataset. Note you will also need the global sampling file to use the auk package for zero filling.

#raw eBird data omits Great Grey Owl & Northern Hawk Owl as sensitive species (https://support.ebird.org/en/support/solutions/articles/48000803210?b_id=1928&_gl=1*xq054u*_ga*ODczMTUyMjcuMTY2OTE0MDI4Ng..*_ga_QR4NVXZ8BM*MTY2OTE0MDI4NS4xLjEuMTY2OTE0MDM3OC4zNS4wLjA.&_ga=2.147122167.150058226.1669140286-87315227.1669140286) and should not be used for modelling these two species.

#wrangling eBird data with the auk package requires installation of AWK on windows computers. Please see #https://cornelllabofornithology.github.io/auk/articles/auk.html.

#eBird data has not been zerofilled because there was no species filtering done and we are assuming that all stationary counts have at least 1 bird observed.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(auk) #eBird wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Set the eBird version ----
v.ebd <- "Jan-2026"

#FILTER DATA###############

#1. Set ebd path----
auk_set_ebd_path(file.path(root, "eBird", v.ebd, paste0("ebd_CA_smp_rel", v.ebd)), overwrite=TRUE)

#2. Define filters----
filters <- auk_ebd(file="ebd_CA_relJan-2026.txt") |>
  auk_protocol("Stationary") |>
  auk_duration(c(1, 10)) |>
  auk_complete() 

#3. Filter data----
#select columns to keep
filtered <- auk_filter(filters, file=file.path(root, "eBird", v.ebd, paste0("03_ebd_filtered_", v.ebd, ".txt")), overwrite=TRUE,
                       keep = c("group identifier", "sampling_event_identifier", "scientific name", "common_name", "observation_count", "latitude", "longitude", "locality_type", "observation_date", "time_observations_started", "observer_id", "duration_minutes"))

