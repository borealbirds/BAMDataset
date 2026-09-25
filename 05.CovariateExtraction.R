# ---
# title: BAM dataset - extract covariates
# author: Elly Knight
# created: September 2, 2026
# ---

#NOTES################################

#PURPOSE: This script extracts point-level covariates from the master BAMSpatialData asset catalogue. It includes annual CanLaD disturbance class, the static NTEMS wildfire dNBR product, and SCANFI biomass, age, and land cover. CanLaD is matched to the survey year; SCANFI uses the nearest available five-year layer.

#Covariates are kept in the R session, checked, and then optionally written as
#one table named covariate in the same BAMDataset DuckDB used to read visits.
#Existing database tables are left unchanged. No manifests, RDS files, or
#standalone CSV outputs are created.

#The current settings run the complete extraction and write the covariate table.
#For a trial run, set test.run to TRUE and write.covariate.table to FALSE.

#PREAMBLE############################

#1. Load packages----
library(DBI) #read the BAMDataset DuckDB
library(duckdb) #connect to the BAMDataset DuckDB
library(BAMCovariates) #build and run covariate extraction jobs

#This script requires BAMCovariates 0.4.0 or later for assets.csv support
if (utils::packageVersion("BAMCovariates") < "0.4.0") {
  stop(
    "BAMCovariates 0.4.0 or later is required. ",
    "Restart R after installing the current package version."
  )
}

#Install a published update with devtools::install_github("borealbirds/BAMCovariates")
#During local development use devtools::install("<path-to-local-BAMCovariates-clone>")

#2. Set root paths----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
spatial.root <- "G:/Shared drives/BAM_SpatialData"

#Locate the BAMSpatialData catalogue independently of the current R working
#directory. This default allows Windows usernames to differ among computers.
#Set BAM_SPATIAL_CATALOG before running if the repository is cloned elsewhere.
spatial.catalog.path <- Sys.getenv(
  "BAM_SPATIAL_CATALOG",
  unset = file.path(
    Sys.getenv("USERPROFILE"),
    "Documents",
    "BAM",
    "Data",
    "BAMSpatialData",
    "assets.csv"
  )
)

#3. Set the BAMDataset versions----
v.wt <- "2026-08-26"
v.ebd <- "Jun-2026"

#4. Choose whether to run a small test or the complete dataset----
#Leave test.run as TRUE until the first extraction has been checked
test.run <- FALSE
test.n.surveys <- 1000

#5. Choose whether to write the covariate table----
#Keep this FALSE while test.run is TRUE. The extraction remains in the covariates object.
write.covariate.table <- TRUE

#Keep this FALSE to protect an existing covariate table
overwrite.covariate.table <- FALSE

#6. Set job size----
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

#2. Locate the master spatial asset catalogue----
if (!file.exists(spatial.catalog.path)) {
  stop(
    "BAMSpatialData assets.csv does not exist: ",
    spatial.catalog.path
  )
}
spatial.catalog.path <- normalizePath(
  spatial.catalog.path,
  winslash = "/",
  mustWork = TRUE
)

#SURVEYS#############################

#1. Connect to the BAMDataset----
database.connection <- dbConnect(
  duckdb(),
  dbdir = database.path,
  read_only = TRUE
)

#2. Read one row per survey----
#Extract the fields required by the package without filtering on the spatial
#extent. The extraction table below applies in_Canada only while reading raster
#values, so the source BAMDataset itself retains all surveys.
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

#COVARIATE EXTRACTION TABLE##############

#1. Select logical raster series from the master asset catalogue----
#The package collapses all physical year files into one extraction row per
#series. Only assets registered as ready by their download scripts are used.
selected.asset.series <- c(
  "canlad_annual_class",
  "ntems_wildfire_dnbr",
  "scanfi_biomass",
  "scanfi_age",
  "scanfi_nfi_landcover"
)

extraction.table <- create_extraction_table(
  assets = spatial.catalog.path,
  asset_series_ids = selected.asset.series
)

#2. Configure every selected series as a point-level extraction----
#create_extraction_table() deliberately returns disabled point-extraction rows.
#Enable them here and give the DuckDB columns concise, stable names.
extraction.table$enabled <- TRUE
extraction.table$statistic <- "value"
extraction.table$buffer_m <- 0

#Enforce the Canadian extent in the extraction specification as well as in the
#DuckDB query. This protects the workflow if its survey query is changed later.
extraction.table$filter_column <- "in_Canada"

output.names <- c(
  canlad_annual_class = "canlad_class",
  ntems_wildfire_dnbr = "wildfire_dnbr_1985_2022",
  scanfi_biomass = "scanfi_biomass",
  scanfi_age = "scanfi_age",
  scanfi_nfi_landcover = "scanfi_landcover"
)
extraction.table$output_name <- unname(
  output.names[extraction.table$asset_series_id]
)

if (
  anyNA(extraction.table$output_name) ||
    !setequal(extraction.table$asset_series_id, selected.asset.series)
) {
  stop("The extraction table does not contain every requested asset series")
}

#3. Inspect the extraction specification----
#CanLaD should be annual, SCANFI nearest, and wildfire dNBR static.
extraction.table

#EXTRACT##############################

#1. Extract all selected point values into the R session----
#Survey coordinates are longitude/latitude (EPSG:4326). buffer_crs is retained for compatibility with future buffer extractions but is not used when buffer_m is zero.
covariates <- extract_covariates(
  surveys = surveys,
  extraction_table = extraction.table,
  raster_root = spatial.root,
  survey_crs = 4326,
  buffer_crs = 3978,
  chunk_size = chunk.size
)

#2. Check the combined covariate table----
nrow(covariates)
head(covariates)
table(covariates$canlad_class, useNA = "ifany")
table(covariates$scanfi_landcover, useNA = "ifany")
summary(
  covariates[c(
    "wildfire_dnbr_1985_2022",
    "scanfi_biomass",
    "scanfi_age"
  )]
)

if (anyDuplicated(covariates$survey_id)) {
  stop("The covariate table contains duplicate survey_id values")
}

#WRITE COVARIATE TABLE################

#1. Protect the test workflow----
#A test extraction should be inspected in the R session and not written to the
#shared BAMDataset.
if (write.covariate.table && test.run) {
  stop("Set test.run to FALSE before writing the covariate table")
}

#2. Reconnect to the existing BAMDataset with write access----
if (write.covariate.table) {
  database.connection <- dbConnect(
    duckdb(),
    dbdir = database.path,
    read_only = FALSE
  )

  tryCatch(
    {
      #3. Protect or replace an existing covariate table----
      table.exists <- dbExistsTable(database.connection, "covariate")
      if (table.exists && !overwrite.covariate.table) {
        stop(
          "The covariate table already exists in: ",
          database.path,
          "\nSet overwrite.covariate.table to TRUE only when replacement ",
          "is intended."
        )
      }

      #4. Write within a transaction----
      #If writing fails, DuckDB rolls the transaction back and retains the
      #previous database state.
      dbWithTransaction(
        database.connection,
        dbWriteTable(
          database.connection,
          name = "covariate",
          value = covariates,
          overwrite = overwrite.covariate.table
        )
      )

      #5. Verify the completed table----
      written.rows <- dbGetQuery(
        database.connection,
        "SELECT COUNT(*) AS n FROM covariate"
      )$n[[1]]

      if (written.rows != nrow(covariates)) {
        stop("The number of rows written to the covariate table is incorrect")
      }

      message("Added covariate table to BAMDataset: ", database.path)
      message("Rows written to covariate table: ", written.rows)
    },
    finally = {
      dbDisconnect(database.connection, shutdown = TRUE)
    }
  )
}
