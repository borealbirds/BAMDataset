# ---
# title: BAM dataset - download WildTrax data
# author: Elly Knight
# created: March 2, 2026
# ---

#NOTES################################

#PURPOSE: This script downloads WildTrax data for the BAM dataset. This script and the subsequent sets of scripts should be run biannually to update the BAM dataset for modelling and use by contributing scientists.

#We use a for loop because it allows for debugging and evaluation of project download status

#ARU main and AI reports are requested together and stored as separate tables.
#Each report is checkpointed by project before the next project is downloaded.
#Rerunning the same snapshot date reuses completed checkpoints, which avoids
#repeating large downloads after an interrupted run.

#There are a handful of projects that are not downloading properly via wildRtrax due to special characters. An issue is open on this. These projects are listed by report in the error.log object. This output should be reviewed and any additional datasets desired should be downloaded manually and placed in the "website downloads (error projects)" directory. Note that many of these datasets are single species call rate datasets that will be of less interest for BAM's purposes.

#Projects that are removed for exclusion are projects for which we have access to the entire organization (e.g., CWS), but where the organization managers would prefer we exclude the projects from the dataset. This is instead of providing access to individual projects within the organization, which can be unwieldy for large organization. Hopefully, an "exclude" button will be available within WildTrax sharing in the future.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(wildrtrax) #to download data from wildtrax
library(readxl) #to read excel files
library(data.table)
library(duckdb)
library(sf)
library(duckspatial)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Choose whether to run a small test----
#A test uses a reproducible random sample of ARU and point-count projects and
#writes to separate checkpoints and a separate DuckDB.
test.run <- TRUE
test.n.projects.per.sensor <- 5
test.seed <- 123

#4. Set the snapshot date----
#Set BAM_WT_SNAPSHOT_DATE to an existing YYYY-MM-DD snapshot to resume a run
#on a later day. If it is unset, a new snapshot is started using today's date.
snapshot.date <- as.Date(Sys.getenv(
  "BAM_WT_SNAPSHOT_DATE",
  unset = as.character(Sys.Date())
))
if (is.na(snapshot.date)) {
  stop("BAM_WT_SNAPSHOT_DATE must use YYYY-MM-DD format")
}

out.dir <- if (test.run) {
  file.path(root, "WildTrax", "test")
} else {
  file.path(root, "WildTrax", as.character(snapshot.date))
}
cache.dir <- file.path(out.dir, "report_cache")
dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)

#AI reports can be much larger than main reports. Allow up to 30 minutes for
#one project download before retrying it.
download.timeout.seconds <- 30 * 60

#5. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#INVENTORY#################

#See printed emails in "Data Assessment/CWS engagement" on the shared drive for further details on the preferences of each CWS region for how to handle data use.

#1. Get project list -----
proj.aru <- wt_get_projects("ARU")
proj.pc <- wt_get_projects("PC")

#2. Get list for exclusion for CWS-PRA ----
ex <- read_excel(file.path(
  "G:/Shared drives/BAM_AvianData/Data Sharing",
  "Dataset Assessment",
  "Exclusion",
  "CWS-PRA_projects_Apr2026.xlsx"
)) |>
  dplyr::filter(Sharing == "Exclude")

#3. Filter project list ----
#remove certain project statuses for CWS-ONT as per preferences
proj <- rbind(proj.aru, proj.pc) |>
  dplyr::filter(
    project_status != "Test Only",
    tasks_completed > 0,
    !project_id %in% ex$project_id,
    !(organization_name == "CWS-ONT" &
      project_status %in% c("Active", "Published = Private"))
  )

#4. Subsample projects for a test run ----
if (test.run) {
  set.seed(test.seed)
  proj <- proj |>
    arrange(project_sensor, project_id) |>
    group_by(project_sensor) |>
    slice_sample(n = test.n.projects.per.sensor) |>
    ungroup()
}

#DOWNLOAD ###############

#1. Set up checkpoint helpers ----
n_proj <- nrow(proj)
error.rows <- list()

checkpoint.path <- function(sensor, project_id, report) {
  file.path(
    cache.dir,
    paste(sensor, project_id, report, sep = "_") |>
      paste0(".parquet")
  )
}

sql.path <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  paste0("'", gsub("'", "''", normalized, fixed = TRUE), "'")
}

parquet.source <- function(paths) {
  if (!length(paths)) {
    stop("At least one Parquet checkpoint is required")
  }
  paste0(
    "read_parquet([",
    paste(vapply(paths, sql.path, character(1)), collapse = ", "),
    "], union_by_name = true)"
  )
}

cache.connection <- dbConnect(duckdb(), dbdir = ":memory:")

write.checkpoint <- function(data, path) {
  if (!is.data.frame(data)) {
    stop("A downloaded report is not a data frame: ", path)
  }

  temporary.path <- tempfile(
    pattern = "checkpoint_",
    tmpdir = dirname(path),
    fileext = ".parquet"
  )
  on.exit(unlink(temporary.path), add = TRUE)

  dbWriteTable(
    cache.connection,
    name = "checkpoint_report",
    value = data,
    overwrite = TRUE
  )
  dbExecute(
    cache.connection,
    paste0(
      "COPY checkpoint_report TO ",
      sql.path(temporary.path),
      " (FORMAT PARQUET, COMPRESSION ZSTD)"
    )
  )

  if (!file.rename(temporary.path, path)) {
    stop("Could not finalize report checkpoint: ", path)
  }
  invisible(path)
}

extract.report <- function(download, report, n_requested) {
  if (is.data.frame(download)) {
    if (n_requested == 1L) {
      return(download)
    }

    #When one of several requested reports is absent, wildrtrax may simplify
    #the remaining report to a data frame. Identify that report from fields
    #that occur only in the AI output so the available table is still saved.
    returned.report <- if (
      any(
        c(
          "confidence",
          "version",
          "is_species_allowed_in_project"
        ) %in%
          names(download)
      )
    ) {
      "ai"
    } else {
      "main"
    }
    if (report != returned.report) {
      stop("The WildTrax download did not contain the ", report, " report")
    }
    return(download)
  }

  if (!is.list(download) || is.null(names(download))) {
    stop("WildTrax returned an unrecognized report object")
  }

  report.pattern <- paste0("_", report, "(_|\\.)")
  matches <- grepl(report.pattern, names(download), ignore.case = TRUE)
  if (!any(matches)) {
    stop("The WildTrax download did not contain the ", report, " report")
  }

  rbindlist(download[matches], fill = TRUE) |>
    as.data.frame()
}

download.with.retries <- function(project_id, sensor, reports, num_times = 5L) {
  last.error <- NULL

  for (attempt in seq_len(num_times)) {
    result <- tryCatch(
      {
        #Long downloads can outlive an authentication token.
        wt_auth()
        wt_download_report(
          project_id = project_id,
          sensor_id = sensor,
          reports = reports,
          max_seconds = download.timeout.seconds
        )
      },
      error = function(e) {
        last.error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(result)) {
      return(list(data = result, attempts = attempt, error = NULL))
    }
  }

  list(data = NULL, attempts = num_times, error = last.error)
}

log.download.error <- function(project_row, report, attempts, message) {
  error.rows[[length(error.rows) + 1L]] <<- project_row |>
    mutate(
      report = report,
      attempts = attempts,
      download_error = message
    )
}

#2. Download reports ----
for (i in seq_len(n_proj)) {
  sensor <- proj$project_sensor[[i]]
  project.id <- proj$project_id[[i]]
  required.reports <- if (sensor == "ARU") c("main", "ai") else "main"
  report.paths <- vapply(
    required.reports,
    function(report) checkpoint.path(sensor, project.id, report),
    character(1)
  )
  missing.reports <- required.reports[!file.exists(report.paths)]

  if (!length(missing.reports)) {
    message(
      "Using cached reports for dataset ",
      proj$project[[i]],
      " : ",
      i,
      " of ",
      n_proj
    )
    next
  }

  downloaded <- download.with.retries(
    project_id = project.id,
    sensor = sensor,
    reports = missing.reports,
    num_times = 5L
  )

  if (is.null(downloaded$data)) {
    for (report in missing.reports) {
      log.download.error(
        project_row = proj[i, ],
        report = report,
        attempts = downloaded$attempts,
        message = downloaded$error
      )
    }
  } else {
    for (report in missing.reports) {
      report.data <- tryCatch(
        extract.report(downloaded$data, report, length(missing.reports)),
        error = function(e) e
      )

      if (inherits(report.data, "error")) {
        log.download.error(
          project_row = proj[i, ],
          report = report,
          attempts = downloaded$attempts,
          message = conditionMessage(report.data)
        )
      } else {
        write.checkpoint(
          report.data,
          checkpoint.path(sensor, project.id, report)
        )
      }
    }
  }

  message(
    "Finished dataset ",
    proj$project[[i]],
    " : ",
    i,
    " of ",
    n_proj
  )
}

#PUT TOGETHER##############

#1. Get manually downloaded reports for errored projects ----
#Note most are published-map only projects that we don't have organization-level access. We don't filter these out earlier because we do have access to some published-map only projects for which we have organization access.
manual.directory <- file.path(
  root,
  "WildTrax",
  "website downloads (error projects)"
)
manual.files <- list.files(
  manual.directory,
  pattern = "\\.csv$",
  full.names = TRUE
)

for (i in seq_along(manual.files)) {
  file.name <- basename(manual.files[[i]])
  file.parts <- str_match(
    file.name,
    regex("^(ARU|PC)_([0-9]+)_(main|ai)(?:_|\\.)", ignore_case = TRUE)
  )

  if (is.na(file.parts[1, 1])) {
    warning("Skipping unrecognized manual report filename: ", file.name)
    next
  }

  sensor <- toupper(file.parts[1, 2])
  project.id <- file.parts[1, 3]
  report <- tolower(file.parts[1, 4])

  #Only incorporate manual reports belonging to this run's project inventory.
  if (
    !any(
      proj$project_sensor == sensor &
        as.character(proj$project_id) == project.id
    )
  ) {
    next
  }

  path <- checkpoint.path(sensor, project.id, report)

  if (file.exists(path)) {
    message("Keeping API checkpoint instead of manual file: ", file.name)
    next
  }

  write.checkpoint(fread(manual.files[[i]]), path)
  message("Checkpointed manual report ", i, " of ", length(manual.files))
}

#2. Find completed report checkpoints ----
aru.main.paths <- list.files(
  cache.dir,
  pattern = "^ARU_[0-9]+_main\\.parquet$",
  full.names = TRUE
)
aru.ai.paths <- list.files(
  cache.dir,
  pattern = "^ARU_[0-9]+_ai\\.parquet$",
  full.names = TRUE
)
pc.main.paths <- list.files(
  cache.dir,
  pattern = "^PC_[0-9]+_main\\.parquet$",
  full.names = TRUE
)

if (!length(aru.main.paths)) {
  stop("No ARU main reports were downloaded")
}
if (!length(aru.ai.paths)) {
  stop("No ARU AI reports were downloaded")
}
if (!length(pc.main.paths)) {
  stop("No point-count main reports were downloaded")
}

#3. Assemble the smaller main reports in memory ----
#The much larger AI reports are copied directly from their Parquet checkpoints
#to the final DuckDB below and are never combined into one in-memory R object.
aru.df <- dbGetQuery(
  cache.connection,
  paste0("SELECT * FROM ", parquet.source(aru.main.paths))
) |>
  as.data.table()
pc.df <- dbGetQuery(
  cache.connection,
  paste0("SELECT * FROM ", parquet.source(pc.main.paths))
) |>
  as.data.table()

aru.df <- aru.df[is.finite(longitude) & is.finite(latitude)]
pc.df <- pc.df[is.finite(longitude) & is.finite(latitude)]

#4. Create the report-level error log ----
if (length(error.rows)) {
  error.log <- bind_rows(error.rows)
} else {
  error.log <- proj[0, ] |>
    mutate(
      report = character(),
      attempts = integer(),
      download_error = character()
    )
}
error.log <- error.log |>
  rowwise() |>
  mutate(
    checkpoint_available = file.exists(checkpoint.path(
      project_sensor,
      project_id,
      report
    ))
  ) |>
  ungroup()

#5. Save date-stamped data and project inventory ----
UTM_CRS <- 4326

aru.sf <- st_as_sf(aru.df, coords = c("longitude", "latitude"), crs = UTM_CRS)
pc.sf <- st_as_sf(pc.df, coords = c("longitude", "latitude"), crs = UTM_CRS)

## write to duckdb ----

db.conn.path <- file.path(
  out.dir,
  if (test.run) {
    "01_wildtrax_raw_test.duckdb"
  } else {
    paste0("01_wildtrax_raw_", snapshot.date, ".duckdb")
  }
)

# make connection
ddbs_write_dataset(aru.sf, db.conn.path, layer = "aru", overwrite = TRUE)

# re-load connection to be able to add more data
db.conn <- dbConnect(duckdb(db.conn.path))
dbExecute(db.conn, "LOAD spatial") # have to do this for whatever reason

# write the other reports and metadata
ddbs_write_table(db.conn, pc.sf, name = "pc", overwrite = TRUE)
dbExecute(
  db.conn,
  paste0(
    "CREATE OR REPLACE TABLE aru_ai AS SELECT * FROM ",
    parquet.source(aru.ai.paths)
  )
)
dbWriteTable(db.conn, name = "project_inventory", proj, overwrite = TRUE)
dbWriteTable(db.conn, name = "error_log", error.log, overwrite = TRUE)

message("WildTrax snapshot written to: ", db.conn.path)
message("Test run: ", test.run)
message("ARU main project reports: ", length(aru.main.paths))
message("ARU AI project reports: ", length(aru.ai.paths))
message("Point-count main project reports: ", length(pc.main.paths))
message("Unresolved report errors: ", sum(!error.log$checkpoint_available))

dbDisconnect(db.conn, shutdown = TRUE)
dbDisconnect(cache.connection, shutdown = TRUE)
