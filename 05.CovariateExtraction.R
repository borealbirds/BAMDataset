# ---
# title: BAM dataset - extract covariates
# author: Elly Knight
# created: September 2, 2026
# updated: September 26, 2026
# ---

#NOTES################################

#PURPOSE: This is the user-facing controller for covariate extraction. The
#same script can run restartable jobs locally, prepare jobs for DRAC, run one
#Slurm array worker, or combine completed jobs and update the BAMDataset
#DuckDB. BAMCovariates contains the GIS implementation; this script handles
#BAMDataset paths, survey input, staging, checks, and database writing.

#Each job processes one survey chunk against every selected covariate. Job RDS
#files and the combined RDS are restart checkpoints, not final data products.
#The final product is the covariate table in the existing BAMDataset DuckDB.
#Keep the checkpoints until that table has been inspected and verified.

#The recipe's filter_column controls the spatial extent of each covariate. The
#current recipe uses in_Canada, so all visits remain in the final table while
#new values outside Canada remain missing.

#WORKFLOW MODES######################

#Choose one mode when running interactively:
#  local         Build jobs, run them locally, combine them, and optionally
#                update DuckDB.
#  drac_prepare  Build the manifest and survey chunks on a DRAC login node.
#                This mode does not run extraction or update DuckDB.
#  combine       Validate completed jobs, combine them, and optionally update
#                DuckDB. Use after the Slurm array has finished.
#  drac_worker   Run one Slurm array job. The submission script selects this
#                mode automatically; users normally do not select it here.
workflow.mode <- "local"

#A command-line mode overrides the value above. For example, the Slurm template
#runs: Rscript 05.CovariateExtraction.R drac_worker /path/job_manifest.csv
script.arguments <- commandArgs(trailingOnly = TRUE)
valid.workflow.modes <- c("local", "drac_prepare", "combine", "drac_worker")
if (length(script.arguments) && script.arguments[[1]] %in% valid.workflow.modes) {
  workflow.mode <- script.arguments[[1]]
} else if (length(script.arguments) && !interactive()) {
  stop(
    "Unknown workflow mode: ", script.arguments[[1]],
    ". Choose one of: ", paste(valid.workflow.modes, collapse = ", ")
  )
}
if (!workflow.mode %in% valid.workflow.modes) {
  stop(
    "workflow.mode must be one of: ",
    paste(valid.workflow.modes, collapse = ", ")
  )
}

#PREAMBLE############################

#1. Load and check BAMCovariates----
library(BAMCovariates) #validate recipes and perform spatial extraction

if (utils::packageVersion("BAMCovariates") < "0.12.0") {
  stop(
    "BAMCovariates 0.12.0 or later is required. ",
    "Restart R after installing the current package version."
  )
}

#Install a published update with:
#devtools::install_github("borealbirds/BAMCovariates")

#During local package development use:
#devtools::install("C:/Users/<username>/Documents/BAM/Data/BAMCovariates")

#DRAC WORKER#########################

#A Slurm worker needs only the saved manifest and BAMCovariates. It does not
#open DuckDB or load the full survey dataset. Keeping this branch first makes
#individual array tasks small and independent.
if (workflow.mode == "drac_worker") {
  manifest.path <- if (length(script.arguments) >= 2L) {
    script.arguments[[2L]]
  } else {
    Sys.getenv("BAM_COVARIATE_MANIFEST", unset = "")
  }
  if (!nzchar(manifest.path) || !file.exists(manifest.path)) {
    stop(
      "The DRAC worker cannot find its manifest. Pass the job_manifest.csv ",
      "path as the second command-line argument or set BAM_COVARIATE_MANIFEST."
    )
  }

  completed.path <- run_job(
    manifest = manifest.path,
    job_id = slurm_job_id(),
    progress = FALSE
  )
  message("DRAC worker completed: ", completed.path)
}

#CONTROLLER SETTINGS#################

#The remaining sections are not needed inside a Slurm worker.
if (workflow.mode != "drac_worker") {

  #1. Load database packages----
  library(DBI) #read and update the BAMDataset DuckDB
  library(duckdb) #connect to the BAMDataset DuckDB

  #2. Set local defaults or override them with environment variables----
  #The G: defaults work on BAM Windows computers. On DRAC, define these three
  #environment variables before running drac_prepare or combine:
  #  BAM_DATASET_ROOT       directory containing the BAMDataset DuckDB
  #  BAM_SPATIAL_ROOT       root of the staged native raster files
  #  BAM_COVARIATE_JOB_ROOT scratch/project directory for restart checkpoints
  root <- Sys.getenv(
    "BAM_DATASET_ROOT",
    unset = "G:/Shared drives/BAM_AvianData/BAMDataset"
  )
  spatial.root <- Sys.getenv(
    "BAM_SPATIAL_ROOT",
    unset = "G:/Shared drives/BAM_SpatialData"
  )
  job.root <- Sys.getenv(
    "BAM_COVARIATE_JOB_ROOT",
    unset = file.path(root, "05_covariate_jobs")
  )

  #Locate this repository independently of the current working directory.
  #Set BAM_DATASET_REPOSITORY if the repository is cloned elsewhere.
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

  #4. Define this extraction stage----
  #Use a short name containing only letters, numbers, underscores, or hyphens.
  #Set stage.outputs to NULL to consider every enabled recipe output. When DRAC
  #storage is limited, list only the output_name values whose rasters are staged.
  stage.name <- "all_covariates"
  stage.outputs <- NULL

  #Examples for storage-limited stages:
  #stage.name <- "forest_structure"
  #stage.outputs <- c("scanfi_biomass", "scanfi_age", "scanfi_landcover")

  #5. Choose whether to run a small test----
  #A test uses only the first test.n.surveys visits, ignores existing
  #covariates, writes to a separate checkpoint directory, and cannot update
  #the shared DuckDB.
  test.run <- FALSE
  test.n.surveys <- 1000L

  #6. Choose whether existing output columns may be replaced----
  #FALSE retains existing DuckDB columns and extracts only requested outputs
  #that are not already present. Set TRUE deliberately to recalculate and
  #replace the outputs selected by stage.outputs.
  replace.existing.outputs <- FALSE

  #7. Choose whether to update the DuckDB covariate table----
  #Leave FALSE until the combined covariates object has been inspected. Only
  #local and combine modes can write. drac_prepare never changes the database.
  write.covariate.table <- FALSE

  #8. Set the maximum surveys per restartable job----
  #Start with 50,000. Reduce this if a worker exceeds its requested memory.
  chunk.size <- if (test.run) test.n.surveys else 50000L

  #9. Choose whether to rebuild the manifest before local/preparation runs----
  #TRUE is recommended. Existing valid result chunks are still reused when the
  #survey chunk, recipe, raster sources, and package version have not changed.
  rebuild.manifest <- TRUE

  #PATHS###############################

  #1. Check the stage name and create its paths----
  if (
    length(stage.name) != 1L || is.na(stage.name) || !nzchar(stage.name) ||
      !grepl("^[A-Za-z0-9_-]+$", stage.name)
  ) {
    stop(
      "stage.name must contain only letters, numbers, underscores, or hyphens"
    )
  }
  run.name <- if (test.run) paste0(stage.name, "_test") else stage.name
  stage.directory <- file.path(job.root, run.name)
  manifest.path <- file.path(stage.directory, "job_manifest.csv")
  combined.path <- file.path(stage.directory, "combined_covariates.rds")

  #An explicit second argument may identify an existing manifest for combine.
  if (
    workflow.mode == "combine" && length(script.arguments) >= 2L &&
      nzchar(script.arguments[[2L]])
  ) {
    manifest.path <- script.arguments[[2L]]
    stage.directory <- dirname(manifest.path)
    combined.path <- file.path(stage.directory, "combined_covariates.rds")
  }

  #2. Locate the BAMDataset DuckDB----
  default.database.path <- file.path(
    root,
    paste0("03_BAMDataset_WT-", v.wt, "_EBd-", v.ebd, ".duckdb")
  )
  database.path <- Sys.getenv(
    "BAM_DATASET_DATABASE",
    unset = default.database.path
  )
  if (!file.exists(database.path)) {
    stop("BAMDataset DuckDB does not exist: ", database.path)
  }

  #3. Locate the committed extraction recipe----
  extraction.table.path <- Sys.getenv(
    "BAM_COVARIATE_RECIPE",
    unset = file.path(repository.root, "05.CovariateExtractionTable.csv")
  )
  if (!file.exists(extraction.table.path)) {
    stop("Covariate extraction table does not exist: ", extraction.table.path)
  }

  #SURVEYS#############################

  #1. Connect to the BAMDataset for reading----
  database.connection <- dbConnect(
    duckdb(),
    dbdir = database.path,
    read_only = TRUE
  )

  #2. Read every survey without applying a spatial filter----
  #Only fields required for extraction are brought into R. ORDER BY makes the
  #survey chunks reproducible. LIMIT applies only to the temporary test run.
  survey.query <- paste(
    "SELECT",
    "  survey_id,",
    "  longitude,",
    "  latitude,",
    "  in_Canada,",
    "  CAST(EXTRACT(YEAR FROM date_time) AS INTEGER) AS survey_year",
    "FROM visit",
    "ORDER BY survey_id"
  )
  if (test.run) {
    survey.query <- paste(survey.query, "LIMIT", test.n.surveys)
  }
  surveys <- dbGetQuery(database.connection, survey.query)

  #3. Read existing covariates for incremental extraction----
  #Test runs deliberately start fresh so they remain small and independent.
  existing.covariates <- NULL
  if (!test.run && dbExistsTable(database.connection, "covariate")) {
    existing.covariates <- dbReadTable(database.connection, "covariate")
  }
  dbDisconnect(database.connection, shutdown = TRUE)

  #4. Check survey identifiers----
  surveys$survey_id <- as.character(surveys$survey_id)
  if (anyNA(surveys$survey_id) || any(!nzchar(surveys$survey_id))) {
    stop("visit contains a blank or missing survey_id")
  }
  if (anyDuplicated(surveys$survey_id)) {
    stop("visit contains duplicate survey_id values")
  }

  #5. Retain usable coordinate rows for raster extraction----
  #Every survey ID is restored when the stage results are joined below.
  usable.coordinates <- !is.na(surveys$longitude) &
    !is.na(surveys$latitude) &
    is.finite(surveys$longitude) &
    is.finite(surveys$latitude)
  extraction.surveys <- surveys[usable.coordinates, , drop = FALSE]
  if (!nrow(extraction.surveys)) {
    stop("No surveys have usable coordinates for covariate extraction")
  }

  #6. Check and subset the existing covariate table----
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
        "First unmatched value: ", unknown.ids[[1L]]
      )
    }
    existing.for.extraction <- existing.covariates[
      existing.covariates$survey_id %in% extraction.surveys$survey_id,
      ,
      drop = FALSE
    ]
  }

  #7. Inspect the survey input----
  message("Surveys in this run: ", nrow(surveys))
  message("Surveys with usable coordinates: ", nrow(extraction.surveys))
  print(table(surveys$in_Canada, useNA = "ifany"))
  print(range(surveys$survey_year, na.rm = TRUE))
  print(utils::head(surveys))

  #COVARIATE EXTRACTION TABLE##########

  #1. Read the committed recipe without requiring every staged raster yet----
  extraction.table <- read_extraction_table(
    extraction_table = extraction.table.path,
    raster_root = spatial.root,
    check_files = FALSE
  )
  enabled.output.names <- unique(
    extraction.table$output_name[extraction.table$enabled]
  )

  #2. Validate this stage's requested output names----
  requested.output.names <- if (is.null(stage.outputs)) {
    enabled.output.names
  } else {
    if (
      !is.character(stage.outputs) || !length(stage.outputs) ||
        anyNA(stage.outputs) || any(!nzchar(stage.outputs)) ||
        anyDuplicated(stage.outputs)
    ) {
      stop("stage.outputs must contain unique, non-blank output_name values")
    }
    unknown.outputs <- setdiff(stage.outputs, enabled.output.names)
    if (length(unknown.outputs)) {
      stop(
        "stage.outputs contains values absent from the enabled recipe: ",
        paste(unknown.outputs, collapse = ", ")
      )
    }
    stage.outputs
  }

  #3. Preview incremental extraction decisions----
  extraction.plan <- plan_covariate_extraction(
    extraction_table = extraction.table,
    existing_covariates = existing.for.extraction,
    survey_crs = 4326,
    buffer_crs = 3978
  )
  print(extraction.plan)

  #4. Select only outputs that this run should calculate----
  #DuckDB cannot retain R metadata attributes. Existing DuckDB columns are
  #therefore skip_unverified and are retained by name unless replacement was
  #explicitly requested above.
  if (replace.existing.outputs) {
    outputs.to.extract <- requested.output.names
  } else {
    outputs.to.extract <- extraction.plan$output_name[
      extraction.plan$output_name %in% requested.output.names &
        extraction.plan$action %in% c("extract_new", "reextract_changed")
    ]
  }
  skipped.outputs <- setdiff(requested.output.names, outputs.to.extract)
  if (length(skipped.outputs)) {
    message(
      "Retaining existing covariate output(s): ",
      paste(skipped.outputs, collapse = ", ")
    )
  }

  #BUILD OR RUN JOBS##################

  stage.covariates <- NULL

  #1. Build the manifest for local or DRAC execution----
  if (workflow.mode %in% c("local", "drac_prepare")) {
    if (!length(outputs.to.extract)) {
      message("No requested covariates require extraction; no jobs were built.")
    } else if (file.exists(manifest.path) && !rebuild.manifest) {
      message("Using existing job manifest: ", manifest.path)
    } else {
      jobs <- build_jobs(
        surveys = extraction.surveys,
        extraction_table = extraction.table,
        output_dir = stage.directory,
        survey_crs = 4326,
        buffer_crs = 3978,
        chunk_size = chunk.size,
        overwrite_manifest = file.exists(manifest.path),
        stage_outputs = outputs.to.extract
      )
      manifest.path <- file.path(stage.directory, "job_manifest.csv")
    }
  }

  #2. Run and combine all jobs on the local computer----
  if (workflow.mode == "local" && length(outputs.to.extract)) {
    run_local(manifest.path, progress = TRUE)
    stage.covariates <- combine_results(manifest.path, combined.path)
  }

  #3. Report the next command after DRAC preparation----
  if (workflow.mode == "drac_prepare" && length(outputs.to.extract)) {
    manifest.rows <- nrow(utils::read.csv(manifest.path, check.names = FALSE))
    message("DRAC preparation is complete.")
    message("Manifest: ", manifest.path)
    message("Submit a Slurm array numbered 1-", manifest.rows, ".")
    message(
      "After every worker succeeds, run this script with workflow.mode = ",
      "'combine'."
    )
  }

  #4. Combine completed local or Slurm jobs----
  if (workflow.mode == "combine") {
    if (!file.exists(manifest.path)) {
      stop("Job manifest does not exist: ", manifest.path)
    }
    stage.covariates <- combine_results(manifest.path, combined.path)
    outputs.to.extract <- setdiff(names(stage.covariates), "survey_id")
  }

  #ASSEMBLE FINAL TABLE###############

  #Only local and combine modes create the in-session final covariate table.
  if (workflow.mode %in% c("local", "combine")) {

    #1. Start with one row for every visit----
    covariates <- data.frame(
      survey_id = surveys$survey_id,
      stringsAsFactors = FALSE
    )

    #2. Retain all existing covariate columns----
    if (!is.null(existing.covariates)) {
      existing.positions <- match(
        covariates$survey_id,
        existing.covariates$survey_id
      )
      for (column in setdiff(names(existing.covariates), "survey_id")) {
        covariates[[column]] <- existing.covariates[[column]][existing.positions]
      }
    }

    #3. Add the completed stage by survey_id----
    if (!is.null(stage.covariates)) {
      stage.covariates$survey_id <- as.character(stage.covariates$survey_id)
      if (anyDuplicated(stage.covariates$survey_id)) {
        stop("The completed stage contains duplicate survey_id values")
      }
      unknown.stage.ids <- setdiff(
        stage.covariates$survey_id,
        covariates$survey_id
      )
      if (length(unknown.stage.ids)) {
        stop(
          "The completed stage contains a survey absent from visit: ",
          unknown.stage.ids[[1L]]
        )
      }

      stage.columns <- setdiff(names(stage.covariates), "survey_id")
      overlapping.columns <- intersect(stage.columns, names(covariates))
      stage.positions <- match(
        stage.covariates$survey_id,
        covariates$survey_id
      )

      #A repeated combine may encounter columns already written to DuckDB. If
      #their values are identical they are safely retained. Different values
      #require explicit replacement permission.
      identical.overlaps <- character()
      for (column in overlapping.columns) {
        existing.values <- covariates[[column]][stage.positions]
        incoming.values <- stage.covariates[[column]]
        if (isTRUE(all.equal(
          existing.values,
          incoming.values,
          check.attributes = FALSE
        ))) {
          identical.overlaps <- c(identical.overlaps, column)
        }
      }
      conflicting.columns <- setdiff(
        overlapping.columns,
        identical.overlaps
      )
      if (length(conflicting.columns) && !replace.existing.outputs) {
        stop(
          "The completed stage conflicts with existing DuckDB column(s): ",
          paste(conflicting.columns, collapse = ", "),
          ". Inspect the difference, then set replace.existing.outputs <- ",
          "TRUE only if replacement is intended."
        )
      }
      if (length(identical.overlaps)) {
        message(
          "Completed values already match DuckDB: ",
          paste(identical.overlaps, collapse = ", ")
        )
      }

      columns.to.add <- setdiff(stage.columns, identical.overlaps)
      completed.positions <- match(
        covariates$survey_id,
        stage.covariates$survey_id
      )
      for (column in columns.to.add) {
        covariates[[column]] <- stage.covariates[[column]][completed.positions]
      }
      attr(covariates, "BAMCovariates_metadata") <- attr(
        stage.covariates,
        "BAMCovariates_metadata"
      )
    }

    #CHECK###############################

    #1. Verify the final table structure----
    if (nrow(covariates) != nrow(surveys)) {
      stop("The covariate table does not contain one row per visit")
    }
    if (anyDuplicated(covariates$survey_id)) {
      stop("The covariate table contains duplicate survey_id values")
    }
    if (!identical(covariates$survey_id, surveys$survey_id)) {
      stop("The covariate table survey order changed unexpectedly")
    }

    #2. Confirm Canada-filtered outputs are missing outside Canada----
    canada.outputs <- unique(extraction.table$output_name[
      extraction.table$enabled &
        !is.na(extraction.table$filter_column) &
        extraction.table$filter_column == "in_Canada"
    ])
    canada.outputs <- intersect(canada.outputs, outputs.to.extract)
    outside.canada <- is.na(surveys$in_Canada) | surveys$in_Canada != 1
    if (
      length(canada.outputs) &&
        any(!is.na(covariates[outside.canada, canada.outputs, drop = FALSE]))
    ) {
      stop("New covariate values were assigned outside the in_Canada extent")
    }

    #3. Inspect the completed object before enabling database writing----
    print(utils::head(covariates))
    print(summary(covariates))
    print(attr(covariates, "BAMCovariates_metadata"))
    message("Final covariate rows: ", nrow(covariates))
    message("Final covariate columns: ", ncol(covariates) - 1L)

    #WRITE COVARIATE TABLE################

    #1. Protect test and no-change workflows----
    if (write.covariate.table && test.run) {
      stop("Set test.run to FALSE before writing the covariate table")
    }
    if (write.covariate.table && is.null(stage.covariates)) {
      message("No completed stage is available; DuckDB was not changed.")
    }

    #2. Replace the combined table inside one database transaction----
    if (write.covariate.table && !is.null(stage.covariates)) {
      database.connection <- dbConnect(
        duckdb(),
        dbdir = database.path,
        read_only = FALSE
      )

      tryCatch(
        {
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

          #3. Re-read database dimensions before reporting success----
          written.rows <- dbGetQuery(
            database.connection,
            "SELECT COUNT(*) AS n FROM covariate"
          )$n[[1L]]
          written.columns <- dbListFields(database.connection, "covariate")
          if (written.rows != nrow(covariates)) {
            stop("The number of rows written to covariate is incorrect")
          }
          if (!identical(written.columns, names(covariates))) {
            stop("The columns written to covariate are incorrect")
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
  }
}
