# ---
# title: BAM dataset - extract covariates
# author: Elly Knight
# created: September 2, 2026
# ---

#NOTES################################

#PURPOSE: This script tests the CovariateExtraction package by extracting the annual CanLaD disturbance class at each Canadian survey point. Survey year is matched to raster year, so (for example) a Canadian survey conducted in 2020 is extracted only from the 2020 CanLaD raster.

#Covariates are kept in the R session, checked, and then written as one table named covariate in a new version of the BAMDataset DuckDB. No manifests, RDS files, or standalone CSV outputs are created.

#The script defaults to 1,000 surveys and does not write a database. Review the test result before enabling the full extraction and database write.

#PREAMBLE############################

#1. Load packages----
library(DBI) #read the BAMDataset DuckDB
library(duckdb) #connect to the BAMDataset DuckDB
library(CovariateExtraction) #build and run covariate extraction jobs

#This script requires CovariateExtraction 0.3.0 or later
if (utils::packageVersion("CovariateExtraction") < "0.3.0") {
  stop(
    "CovariateExtraction 0.3.0 or later is required. ",
    "Restart R after installing the current package version."
  )
}

#Install a published update with devtools::install_github("borealbirds/CovariateExtraction")
#During local development use devtools::install("C:/Users/elly/Documents/BAM/Data/CovariateExtraction")

#2. Set root paths----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
spatial.root <- "G:/Shared drives/BAM_SpatialData"

#3. Set the BAMDataset versions----
v.wt <- "2026-08-26"
v.ebd <- "Jun-2026"

#4. Set the CanLaD version----
v.canlad <- "v1.1_20260508"
canlad.folder <- "TS_FinalFilter_20260508"

#5. Choose whether to run a small test or the complete dataset----
#Leave test.run as TRUE until the first extraction has been checked
test.run <- FALSE
test.n.surveys <- 1000

#6. Choose whether to write the new DuckDB----
#Keep this FALSE while test.run is TRUE. The extraction remains in the covariates object.
write.database <- FALSE

#Keep this FALSE to protect an existing 05_ database
overwrite.database <- FALSE

#7. Set job size----
#The full workflow can start with 50,000 surveys per job; reduce this if a job uses too much memory
chunk.size <- if (test.run) test.n.surveys else 50000

#PATHS###############################

#1. Locate the BAMDataset DuckDB----
database.path <- file.path(
  root,
  paste0("03_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, ".duckdb")
)

if (!file.exists(database.path)) {
  stop("BAMDataset DuckDB does not exist: ", database.path)
}

#2. Set the new BAMDataset path----
#The source database is copied before the covariate table is added, leaving the 03_ database unchanged
output.database.path <- file.path(
  root,
  paste0("05_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, ".duckdb")
)

#SURVEYS#############################

#1. Connect to the BAMDataset----
database.connection <- dbConnect(
  duckdb(),
  dbdir = database.path,
  read_only = TRUE
)

#2. Read one row per survey----
#Extract only the fields required by the package. The date filter excludes surveys outside the years supplied by this CanLaD release.
survey.query <- paste(
  "SELECT",
  "  survey_id,",
  "  longitude,",
  "  latitude,",
  "  in_Canada,",
  "  CAST(EXTRACT(YEAR FROM date_time) AS INTEGER) AS survey_year",
  "FROM visit",
  "WHERE date_time >= '1984-01-01'",
  "  AND date_time < '2026-01-01'",
  "  AND longitude IS NOT NULL",
  "  AND latitude IS NOT NULL",
  "ORDER BY survey_id"
)

#Limit the initial trial without changing the query used by the full run
if (test.run) {
  survey.query <- paste(survey.query, "LIMIT", test.n.surveys)
}

surveys <- dbGetQuery(database.connection, survey.query)

#3. Disconnect from the BAMDataset----
dbDisconnect(database.connection, shutdown = TRUE)

#4. Check the survey input----
#The package will perform more detailed validation before extraction
nrow(surveys)
range(surveys$survey_year)
head(surveys)

#CANLAD EXTRACTION TABLE##############

#NOTE: This will get replaced by a table that's built outside of R when we expand to more than just one test dataset

#1. Describe the complete annual raster series with one row----
#Relative paths are resolved against spatial.root by extract_covariates()
extraction.table <- data.frame(
  enabled = TRUE,
  category = "Disturbance",
  covariate_id = "canlad_annual_class",
  raster_path = file.path(
    "raw",
    "canlad",
    v.canlad,
    canlad.folder,
    "canlad_annual_{year}_v1_1_20260508.tif"
  ),
  output_name = "canlad_class",
  statistic = "value",
  buffer_m = 0,
  band = 1,
  year_offset = 0,
  temporal_match = "annual",
  filter_column = "in_Canada",
  na_action = "keep",
  fill_value = NA_real_,
  value_type = "categorical",
  units = "CanLaD class",
  notes = "Annual CanLaD disturbance class extracted at Canadian survey points",
  stringsAsFactors = FALSE
)

#2. Inspect the extraction specification----
extraction.table

#EXTRACT##############################

#1. Extract CanLaD values into the R session----
#Survey coordinates are longitude/latitude (EPSG:4326). buffer_crs is retained for compatibility with future buffer extractions but is not used when buffer_m is zero.
covariates <- extract_covariates(
  surveys = surveys,
  extraction_table = extraction.table,
  raster_root = spatial.root,
  survey_crs = 4326,
  buffer_crs = 3978,
  chunk_size = chunk.size
)

#2. Check the single covariate table----
nrow(covariates)
head(covariates)
table(covariates$canlad_class, useNA = "ifany")

if (anyDuplicated(covariates$survey_id)) {
  stop("The covariate table contains duplicate survey_id values")
}

#WRITE BAM DATASET####################

#1. Protect the test workflow----
#A test extraction should be inspected in the R session and never written into a versioned BAMDataset
if (write.database && test.run) {
  stop("Set test.run to FALSE before writing the new BAMDataset")
}

#2. Copy the existing BAMDataset----
#This creates a complete 05_ database while preserving the existing 03_ database
if (write.database) {
  if (file.exists(output.database.path) && !overwrite.database) {
    stop(
      "Output database already exists: ",
      output.database.path,
      "\nSet overwrite.database to TRUE only when replacement is intended."
    )
  }

  database.copied <- file.copy(
    from = database.path,
    to = output.database.path,
    overwrite = overwrite.database
  )

  if (!database.copied) {
    stop("Could not copy the BAMDataset to: ", output.database.path)
  }

  #3. Add the covariate table to the new database----
  output.database.connection <- dbConnect(
    duckdb(),
    dbdir = output.database.path
  )

  dbWriteTable(
    output.database.connection,
    name = "covariate",
    value = covariates,
    overwrite = TRUE
  )

  #4. Verify and close the new database----
  written.rows <- dbGetQuery(
    output.database.connection,
    "SELECT COUNT(*) AS n FROM covariate"
  )$n[[1]]

  dbDisconnect(output.database.connection, shutdown = TRUE)

  if (written.rows != nrow(covariates)) {
    stop("The number of rows written to the covariate table is incorrect")
  }

  message("Created BAMDataset: ", output.database.path)
  message("Rows written to covariate table: ", written.rows)
}
