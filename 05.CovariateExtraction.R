# ---
# title: BAM dataset - extract covariates
# author: Elly Knight
# created: September 2, 2026
# updated: September 26, 2026
# ---

#NOTES################################

#PURPOSE: This script extracts point-level covariates using the committed
#05.CovariateExtractionTable.csv recipe. The initial recipe includes annual
#CanLaD disturbance class, static NTEMS wildfire dNBR, and nearest-year SCANFI
#biomass, age, and land cover.

#All visits are retained in the final covariate table. Raster extraction is
#limited by the recipe's in_Canada filter. Visits outside Canada or without
#usable coordinates therefore remain in the table with missing new covariates.

#If a covariate table already exists, its columns are supplied to
#BAMCovariates. Existing columns are retained and only newly requested recipe
#outputs are extracted. DuckDB does not preserve R attributes, so columns read
#back from DuckDB are reported as unverified rather than provenance-verified.

#No RDS, CSV, manifest, or standalone covariate output is created. The result
#remains in the covariates object and can be written back to the same DuckDB.

#PREAMBLE############################

#1. Load packages----
library(DBI) #read and update the BAMDataset DuckDB
library(duckdb) #connect to the BAMDataset DuckDB
library(BAMCovariates) #validate recipes and extract covariates

#This workflow uses the validation and incremental-extraction behavior added
#in BAMCovariates 0.10.1.
if (utils::packageVersion("BAMCovariates") < "0.10.1") {
  stop(
    "BAMCovariates 0.10.1 or later is required. ",
    "Restart R after installing the current package version."
  )
}

#Install a published update with:
#devtools::install_github("borealbirds/BAMCovariates")

#During local package development use:
#devtools::install("C:/Users/<username>/Documents/BAM/Data/BAMCovariates")

#2. Set root paths----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
spatial.root <- "G:/Shared drives/BAM_SpatialData"

#Locate this repository independently of the current R working directory.
#Set BAM_DATASET_REPOSITORY before running if the repository is cloned elsewhere.
repository.root <- Sys.getenv(
  "BAM_DATASET_REPOSITORY",
  unset = file.path(
    Sys.getenv("USERPROFILE"),
    "Documents",
    "BAM",
    "Data",
    "BAMDataset"
  )
)

#3. Set the BAMDataset versions----
v.wt <- "2026-08-26"
v.ebd <- "Jun-2026"

#4. Choose whether to run a small test or the complete dataset----
#A test run reads only the first test.n.surveys visits, ignores any existing
#covariate table, and cannot be written to the shared database.
test.run <- FALSE
test.n.surveys <- 1000

#5. Choose whether to update the DuckDB covariate table----
#Leave this FALSE until the covariates object has been inspected. When TRUE,
#the script creates the covariate table or replaces it with the combined table.
write.covariate.table <- TRUE

#6. Set extraction chunk size----
#Start with 50,000 surveys per chunk. Reduce this if extraction uses too much
#memory. This is not the number of rows written to DuckDB at once.
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

#2. Locate the committed extraction recipe----
extraction.table.path <- file.path(
  repository.root,
  "05.CovariateExtractionTable.csv"
)

if (!file.exists(extraction.table.path)) {
  stop("Covariate extraction table does not exist: ", extraction.table.path)
}

#SURVEYS#############################

#1. Connect to the BAMDataset----
database.connection <- dbConnect(
  duckdb(),
  dbdir = database.path,
  read_only = TRUE
)

#2. Read every survey without applying date or spatial filters----
#Only the columns required for extraction are brought into R. In a test run,
#LIMIT restricts this temporary trial; the complete run has no WHERE or LIMIT.
survey.query <- paste(
  "SELECT",
  "  survey_id,",
  "  longitude,",
  "  latitude,",
  "  in_Canada,",
  "  CAST(EXTRACT(YEAR FROM date_time) AS INTEGER) AS survey_year",
  "FROM visit"
)

if (test.run) {
  survey.query <- paste(
    survey.query,
    "ORDER BY survey_id LIMIT",
    test.n.surveys
  )
}

surveys <- dbGetQuery(database.connection, survey.query)

#3. Read existing covariates for incremental extraction----
#Test runs deliberately start fresh so they remain small and easy to inspect.
existing.covariates <- NULL
if (!test.run && dbExistsTable(database.connection, "covariate")) {
  existing.covariates <- dbReadTable(database.connection, "covariate")
}

#4. Disconnect from the read-only database connection----
dbDisconnect(database.connection, shutdown = TRUE)

#5. Check the survey identifiers----
surveys$survey_id <- as.character(surveys$survey_id)
if (anyNA(surveys$survey_id) || any(!nzchar(surveys$survey_id))) {
  stop("visit contains a blank or missing survey_id")
}
if (anyDuplicated(surveys$survey_id)) {
  stop("visit contains duplicate survey_id values")
}

#6. Retain only usable coordinate rows for raster extraction----
#All survey IDs are restored below. This subset prevents missing or infinite
#coordinates from entering spatial operations.
usable.coordinates <- !is.na(surveys$longitude) &
  !is.na(surveys$latitude) &
  is.finite(surveys$longitude) &
  is.finite(surveys$latitude)

extraction.surveys <- surveys[usable.coordinates, , drop = FALSE]
if (!nrow(extraction.surveys)) {
  stop("No surveys have usable coordinates for covariate extraction")
}

#7. Check and subset an existing covariate table----
existing.for.extraction <- NULL
if (!is.null(existing.covariates)) {
  existing.covariates$survey_id <- as.character(existing.covariates$survey_id)
  if (anyDuplicated(existing.covariates$survey_id)) {
    stop("The existing covariate table contains duplicate survey_id values")
  }
  unknown.ids <- setdiff(existing.covariates$survey_id, surveys$survey_id)
  if (length(unknown.ids)) {
    stop(
      "The existing covariate table contains survey IDs absent from visit. ",
      "First unmatched value: ",
      unknown.ids[[1]]
    )
  }
  existing.for.extraction <- existing.covariates[
    existing.covariates$survey_id %in% extraction.surveys$survey_id,
    ,
    drop = FALSE
  ]
}

#8. Inspect the survey input----
nrow(surveys)
sum(usable.coordinates)
table(surveys$in_Canada, useNA = "ifany")
range(surveys$survey_year, na.rm = TRUE)
head(surveys)

#COVARIATE EXTRACTION TABLE##########

#1. Read and validate the committed recipe----
#Edit 05.CovariateExtractionTable.csv to add, disable, or document covariates.
#BAMSpatialData/assets.csv is the starting point when adding a new data source,
#but the analysis recipe is not regenerated every time this script runs.
extraction.table <- read_extraction_table(
  extraction_table = extraction.table.path,
  raster_root = spatial.root,
  check_files = FALSE
)

#2. Inspect the extraction specification----
extraction.table

#3. Preview what BAMCovariates will extract or retain----
extraction.plan <- plan_covariate_extraction(
  extraction_table = extraction.table,
  existing_covariates = existing.for.extraction,
  survey_crs = 4326,
  buffer_crs = 3978
)
extraction.plan

#EXTRACT##############################

#1. Extract requested point values into the R session----
#The package checks coordinates, temporal recipes, raster headers, requested
#bands, and temporal coverage before or during extraction.
extracted.covariates <- extract_covariates(
  surveys = extraction.surveys,
  extraction_table = extraction.table,
  survey_crs = 4326,
  buffer_crs = 3978,
  chunk_size = chunk.size,
  existing_covariates = existing.for.extraction
)

#2. Start the final table with every survey_id----
covariates <- data.frame(
  survey_id = surveys$survey_id,
  stringsAsFactors = FALSE
)

#3. Retain existing values for every survey----
#This includes existing values for visits that cannot enter spatial extraction.
if (!is.null(existing.covariates)) {
  existing.positions <- match(
    covariates$survey_id,
    existing.covariates$survey_id
  )
  for (column in setdiff(names(existing.covariates), "survey_id")) {
    covariates[[column]] <- existing.covariates[[column]][existing.positions]
  }
}

#4. Add new or updated extraction values----
extracted.positions <- match(
  covariates$survey_id,
  extracted.covariates$survey_id
)
for (column in setdiff(names(extracted.covariates), "survey_id")) {
  covariates[[column]] <- extracted.covariates[[column]][extracted.positions]
}

#Keep BAMCovariates metadata attached in the R session. Standard DuckDB tables
#cannot store this R attribute, so it will not persist after database writing.
attr(covariates, "BAMCovariates_metadata") <- attr(
  extracted.covariates,
  "BAMCovariates_metadata"
)

#CHECK###############################

#1. Check the combined covariate table----
nrow(covariates)
head(covariates)
summary(covariates)
attr(covariates, "BAMCovariates_metadata")

if (nrow(covariates) != nrow(surveys)) {
  stop("The covariate table does not contain one row per visit")
}
if (anyDuplicated(covariates$survey_id)) {
  stop("The covariate table contains duplicate survey_id values")
}
if (!identical(covariates$survey_id, surveys$survey_id)) {
  stop("The covariate table survey order changed unexpectedly")
}

#2. Confirm that newly extracted values are limited to Canada----
new.output.names <- extraction.plan$output_name[
  extraction.plan$action %in% c("extract_new", "reextract_changed")
]
outside.canada <- is.na(surveys$in_Canada) | surveys$in_Canada != 1
if (
  length(new.output.names) &&
    any(!is.na(covariates[outside.canada, new.output.names, drop = FALSE]))
) {
  stop("New covariate values were assigned outside the in_Canada extent")
}

#WRITE COVARIATE TABLE################

#1. Protect the test workflow----
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
      #3. Create or replace the combined covariate table in a transaction----
      #Existing columns have already been retained in covariates, so replacing
      #the database table does not discard them.
      table.exists <- dbExistsTable(database.connection, "covariate")
      dbWithTransaction(
        database.connection,
        dbWriteTable(
          database.connection,
          name = "covariate",
          value = covariates,
          overwrite = table.exists
        )
      )

      #4. Verify the completed database table----
      written.rows <- dbGetQuery(
        database.connection,
        "SELECT COUNT(*) AS n FROM covariate"
      )$n[[1]]
      written.columns <- dbListFields(database.connection, "covariate")

      if (written.rows != nrow(covariates)) {
        stop("The number of rows written to the covariate table is incorrect")
      }
      if (!identical(written.columns, names(covariates))) {
        stop("The columns written to the covariate table are incorrect")
      }

      message("Updated covariate table in BAMDataset: ", database.path)
      message("Rows written: ", written.rows)
      message("Covariate columns: ", ncol(covariates) - 1L)
    },
    finally = {
      dbDisconnect(database.connection, shutdown = TRUE)
    }
  )
}
