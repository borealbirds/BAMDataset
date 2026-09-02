# ---
# title: BAM dataset - extract covariates
# author: Elly Knight
# created: September 2, 2026
# ---

#NOTES################################

#PURPOSE: This script tests the CovariateExtraction package by extracting the annual CanLaD disturbance class at each Canadian survey point. Survey year is matched to raster year, so (for example) a Canadian survey conducted in 2020 is extracted only from the 2020 CanLaD raster.

#This first version creates a standalone covariate output keyed by survey_id. It does not modify the BAMDataset DuckDB or join covariates to the visit table.

#The script defaults to 1,000 surveys and one extraction job. Review the test result before enabling the full run.

#PREAMBLE############################

#1. Load packages----
library(DBI) #read the BAMDataset DuckDB
library(duckdb) #connect to the BAMDataset DuckDB
library(devtools) #load the development version of CovariateExtraction

#Load the package directly from the neighbouring repository while it is under development
devtools::load_all(file.path("..", "CovariateExtraction"))

#2. Set root paths----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"
spatial.root <- "G:/Shared drives/BAM_SpatialData"

#3. Set the BAMDataset versions----
v.wt <- "2026-08-26"
v.ebd <- "Jun-2026"

#4. Set the CanLaD version----
v.canlad <- "v1.1_20260508"
canlad.folder <- "TS_FinalFilter_20260508"
canlad.years <- 1984:2025

#5. Choose whether to run a small test or the complete dataset----
#Leave test.run as TRUE until the first extraction has been checked
test.run <- TRUE
test.n.surveys <- 1000

#6. Choose which stages to run----
#The safe defaults create a manifest and run only its first job
create.manifest <- TRUE
overwrite.manifest <- FALSE
run.first.job <- TRUE
run.all.jobs <- FALSE
combine.outputs <- FALSE

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

#2. Set the output directory----
#Keep the test and full-run manifests separate because they refer to different survey inputs
run.type <- if (test.run) "test" else "full"
output.directory <- file.path(
  root,
  "Covariates",
  "CanLaD",
  v.canlad,
  run.type
)
manifest.path <- file.path(output.directory, "job_manifest.csv")

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
#The package will perform more detailed validation when the manifest is built
nrow(surveys)
range(surveys$survey_year)
head(surveys)

#CANLAD EXTRACTION TABLE##############

#1. Build one extraction row for each annual raster----
#Relative paths are resolved against spatial.root by build_jobs()
extraction.table <- data.frame(
  covariate_id = paste0("canlad_class_", canlad.years),
  raster_path = file.path(
    "raw",
    "canlad",
    v.canlad,
    canlad.folder,
    paste0(
      "canlad_annual_",
      canlad.years,
      "_v1_1_20260508.tif"
    )
  ),
  output_name = "canlad_class",
  statistic = "value",
  buffer_m = 0,
  band = 1,
  raster_year = canlad.years,
  filter_column = "in_Canada",
  value_type = "categorical",
  units = "CanLaD class",
  notes = "Annual CanLaD disturbance class extracted at Canadian survey points",
  stringsAsFactors = FALSE
)

#2. Inspect the extraction specification----
head(extraction.table)
tail(extraction.table)

#BUILD JOBS###########################

#1. Create or read the job manifest----
#Survey coordinates are longitude/latitude (EPSG:4326). buffer_crs is retained for compatibility with future buffer extractions but is not used when buffer_m is zero.
if (create.manifest) {
  jobs <- build_jobs(
    surveys = surveys,
    extraction_table = extraction.table,
    output_dir = output.directory,
    raster_root = spatial.root,
    survey_crs = 4326,
    buffer_crs = 3978,
    chunk_size = chunk.size,
    overwrite_manifest = overwrite.manifest
  )
} else {
  jobs <- read.csv(manifest.path, stringsAsFactors = FALSE)
}

#2. Inspect the jobs before extraction----
nrow(jobs)
head(jobs)

#EXTRACT##############################

#1. Run the first job----
#Completed jobs are skipped, so leaving this enabled is safe when the script is rerun
if (run.first.job) {
  run_job(jobs, job_id = jobs$job_id[[1]])
}

#2. Run every job locally----
#Enable only after the first chunk has been inspected
if (run.all.jobs) {
  run_local(jobs)
}

#COMBINE##############################

#1. Combine all completed chunks----
#This requires every job in the manifest to be complete
if (combine.outputs) {
  covariates <- combine_results(
    jobs,
    file.path(output.directory, "canlad_point_covariates.csv")
  )

  nrow(covariates)
  table(covariates$canlad_class, useNA = "ifany")
}
