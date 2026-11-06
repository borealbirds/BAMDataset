library(tidyverse)
library(RTMB)
library(sf)
library(duckdb)
library(duckspatial)
library(terra)
library(mvtnorm)
library(dtplyr)
library(foreach)
library(doParallel)

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

# Folder for per-species offset files (one parquet file per species), stamped with the data versions so
# files from different dataset versions are never mixed. These are assembled into a single "offsets"
# table at the end of the script.
offset_dir = file.path(ifelse(FIR, "/scratch/pt1/bam/offsets_parquet", "offsets_parquet"),
                       paste0("WT-", v.wt, "_EBd-", v.ebd, "_PRISM-", v.prism))
dir.create(offset_dir, recursive = TRUE, showWarnings = FALSE)

# If TRUE, species that already have a file in offset_dir are not recalculated (useful for resuming a
# run that stopped partway). Leave FALSE if anything upstream (data, covariates, fits, code) has changed.
SKIP_EXISTING_OFFSETS = FALSE

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

# Set up parallelization
ncores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")) * as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
if (is.na(ncores)) ncores = 1
registerDoParallel(cores = ncores)

# each species' offsets are written to their own parquet file; the loop returns the file paths
offset_files = foreach(sp = all_species[1:10]) %dopar% {
  
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

    this_offset_file = file.path(offset_dir, paste0("offsets_", sp, ".parquet"))
    if (SKIP_EXISTING_OFFSETS && file.exists(this_offset_file)) {
      message("Offsets already exist for ", sp, "; skipping")
      return(this_offset_file)
    }

    ## Read in the appropriate QPAD fit ----

    this_fit_filepath = search_for_file(qpad_fit_dir, paste0(sp, ".rds"))[1]
    if (is.na(this_fit_filepath)) {
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
      
      prism_conversion_rates = log(pi * visit_cov$distance_km ^ 2 * 100) # area of PRISM plots in hectares

      message("Calculating offsets and metrics...")

      # one row per survey (in visit_cov order); linear predictors for all N_BOOT draws are built inside
      all_par_offsets = get_offset_boot_mat(dat = visit_cov,
                                            r_pars = r_pars,
                                            include_road = !is.na(include_road) && include_road,
                                            tod_to_keep = tod_to_keep,
                                            all_par_names = ALL_POSSIBLE_PAR_NAMES,
                                            chunk_size = 100) %>%
        # correct for PRISM
        mutate(mean_log_offset = ifelse(visit_cov$sensor == "PRISM", prism_conversion_rates, mean_log_offset),
               sd_log_offset = ifelse(visit_cov$sensor == "PRISM", 0, sd_log_offset))

      message("Writing offsets to ", this_offset_file)

      # long format: one row per survey for this species
      write_parquet_duckdb(data.frame(survey_id = visit_cov$survey_id,
                                      species_code = sp,
                                      offset_mean = all_par_offsets$mean_log_offset,
                                      offset_sd = all_par_offsets$sd_log_offset),
                           this_offset_file)

      return(this_offset_file)
    }
  }

  NULL

}

offset_files = unlist(offset_files)
offset_species = str_match(basename(offset_files), "^offsets_(.*)\\.parquet$")[, 2]

# Assemble the per-species files into one wide "offsets" table in the BAM duckdb ----
# Columns are offset_mean_<sp>, offset_sd_<sp> for each species, then survey_id (one row per survey, in
# visit_cov order). Every file has the same surveys in the same order, so DuckDB stitches them together
# side by side (POSITIONAL JOIN) straight from the parquet files in one streaming pass; the full table
# never has to fit in R's memory. An "ids_match" column records whether every file's survey_id agrees
# with the first file's on that row; it is checked and then dropped below.

offset_alias = paste0("f", seq_along(offset_files))
offset_select = c(paste0(offset_alias, ".offset_mean AS ", dbQuoteIdentifier(ANSI(), paste0("offset_mean_", offset_species)), ",\n       ",
                         offset_alias, ".offset_sd AS ", dbQuoteIdentifier(ANSI(), paste0("offset_sd_", offset_species))),
                  "f1.survey_id",
                  paste0("coalesce(", paste(c("true", paste0(offset_alias[-1], ".survey_id = f1.survey_id")), collapse = " AND "), ", false) AS ids_match"))
offset_from = paste0("read_parquet(", dbQuoteString(ANSI(), normalizePath(offset_files, winslash = "/")), ") ", offset_alias,
                     collapse = "\nPOSITIONAL JOIN ")

offsets_query = paste0("CREATE OR REPLACE TABLE offsets_new AS\n",
                       "SELECT ", paste(offset_select, collapse = ",\n       "), "\n",
                       "FROM ", offset_from)

# Write the offsets table to the duckdb database ----

# if there's already a lock file, just run this again and again and wait. I don't think we need the extra Sys.sleep() anymore.
while(file.exists(duckdb_lock_file)) Sys.sleep(0.5)

# once the file no longer exists, write it
write_csv(data.frame(0), duckdb_lock_file)

db_bam_conn = dbConnect(duckdb(db_bam_conn_path))
dbExecute(db_bam_conn, "LOAD spatial") # have to do this for whatever reason

# if this query needs more memory than is available, DuckDB spills to disk; keep that on fast local storage
dbExecute(db_bam_conn, paste0("SET temp_directory = ", dbQuoteString(db_bam_conn, normalizePath(file.path(offset_dir, "duckdb_tmp"), winslash = "/", mustWork = FALSE))))
if (length(offset_files) > 0) {
  dbExecute(db_bam_conn, offsets_query)
  n_mismatch = dbGetQuery(db_bam_conn, "SELECT count(*) AS n FROM offsets_new WHERE NOT ids_match")$n
  if (n_mismatch > 0) {
    # leave any existing offsets table untouched
    dbExecute(db_bam_conn, "DROP TABLE offsets_new")
    dbDisconnect(db_bam_conn)
    file.remove(duckdb_lock_file)
    stop(n_mismatch, " rows had mismatched survey_ids across the per-species offset files (were some files made ",
         "from a different visit_cov, e.g. with SKIP_EXISTING_OFFSETS = TRUE?). The offsets table was not changed.")
  }
  dbExecute(db_bam_conn, "ALTER TABLE offsets_new DROP COLUMN ids_match")
  dbExecute(db_bam_conn, "DROP TABLE IF EXISTS offsets")
  dbExecute(db_bam_conn, "ALTER TABLE offsets_new RENAME TO offsets")
} else {
  message("No offset files were produced; the offsets table was not written")
}
dbDisconnect(db_bam_conn)

# remove the lock so other processes can run
file.remove(duckdb_lock_file)
