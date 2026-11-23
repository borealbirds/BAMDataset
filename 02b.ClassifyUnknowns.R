# ---
# title: BAM dataset - classify unknown ARU detections
# author: Elly Knight
# created: November 20, 2026
# ---

#NOTES################################

#PURPOSE: This script uses classifier results to build a decision table for
#reclassifying unknown bird annotations in acoustic recordings. It does not
#change the raw WildTrax data. Accepted decisions will be incorporated when
#the BAM dataset is harmonized in script 03.

#Unknown categories and candidate species are formatted from the current
#WildTrax species list. Candidate classifier detections must occur within one
#second of the human annotation, be enabled in the project, exceed the minimum
#score, and exceed the next-highest score by the minimum score difference.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(wildrtrax) #species list
library(DBI)
library(duckdb)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_AvianData/BAMDataset"

#3. Choose whether to run a small test----
#A test uses the WildTrax test database and writes to a separate test DuckDB.
test.run <- FALSE

#4. Set the WildTrax version----
v.wt <- "2026-11-20"

#5. Set classifier decision values----
classifier.version <- "Hawkears v1.0.8"
time.buffer.seconds <- 1
minimum.score <- 0.7
minimum.score.difference <- 0.2

#6. Login to WildTrax----
source("WTlogin.R")
wt_auth()

#HELPER FUNCTIONS####################

sql.path <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  sql.string(normalized)
}

sql.string <- function(value) {
  paste0("'", gsub("'", "''", value, fixed = TRUE), "'")
}

#The scientific names attached to unknown categories in WildTrax are used
#where possible to restrict a category to a genus or species pair. These
#corrections address spelling variants in the current WildTrax species object.
scientific.candidates <- function(scientific.name, spp) {
  scientific.name <- scientific.name |>
    str_squish() |>
    str_replace_all(c(
      "Dendrapagus" = "Dendragapus",
      "Aechmorphus" = "Aechmophorus",
      "Melopsiza" = "Melospiza",
      "towsendii" = "townsendi",
      "townsendii" = "townsendi",
      "Falcon sp\\." = "Falco"
    ))

  if (
    is.na(scientific.name) ||
      scientific.name == "" ||
      scientific.name == "Waterfowl"
  ) {
    return(character())
  }

  parts <- str_split(scientific.name, "\\s*/\\s*")[[1]] |>
    str_remove("\\s+sp\\.$") |>
    str_squish()
  first.genus <- word(parts[[1]], 1)

  targets <- map_chr(parts, function(part) {
    if (str_count(part, "\\s+") == 0 && length(parts) > 1) {
      if (str_detect(part, "^[a-z]")) {
        return(paste(first.genus, part))
      }
    }
    part
  })

  exact.matches <- spp |>
    dplyr::filter(scientific_name_key %in% str_to_upper(targets)) |>
    pull(species_code)

  genus.targets <- targets[str_count(targets, "\\s+") == 0]
  genus.matches <- spp |>
    dplyr::filter(scientific_genus %in% genus.targets) |>
    pull(species_code)

  unique(c(exact.matches, genus.matches))
}

candidate.species <- function(
  annotation.code,
  normalized.code,
  unknown.order,
  unknown.family,
  unknown.scientific.name,
  spp
) {
  #UNQK is treated as UNDU. Geese, brant, and swans are not duck candidates.
  if (normalized.code == "UNDU") {
    return(
      spp |>
        dplyr::filter(
          species_family == "Anatidae",
          !str_detect(
            species_common_name,
            regex(
              "goose|brant|swan",
              ignore_case = TRUE
            )
          )
        ) |>
        pull(species_code)
    )
  }

  #An unidentified trill is restricted to warblers and sparrows.
  if (annotation.code == "UNTR") {
    return(
      spp |>
        dplyr::filter(
          species_family %in%
            c(
              "Parulidae",
              "Passerellidae"
            )
        ) |>
        pull(species_code)
    )
  }

  #These broader groups do not contain enough taxonomy in the WildTrax row to
  #derive their candidates without a small interpretation rule.
  if (annotation.code == "UNBP") {
    return(
      spp |>
        dplyr::filter(
          species_order %in%
            c(
              "Accipitriformes",
              "Falconiformes"
            )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "URPT") {
    return(
      spp |>
        dplyr::filter(
          species_order %in%
            c(
              "Accipitriformes",
              "Falconiformes",
              "Strigiformes"
            )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNIWBI") {
    return(
      spp |>
        dplyr::filter(
          species_family %in%
            c(
              "Ardeidae",
              "Threskiornithidae",
              "Threskiornithinae",
              "Ciconiidae"
            )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNSH") {
    return(
      spp |>
        dplyr::filter(
          species_family %in%
            c(
              "Scolopacidae",
              "Charadriidae",
              "Recurvirostridae",
              "Haematopodidae",
              "Jacanidae"
            )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UPEE") {
    return(
      spp |>
        # species that could (liberally) be considered "peeps"
        dplyr::filter(species_code %in% c("LESA", "SESA", "WESA", "WRSA", "BASA", "LIST", "BBIS", "SAND", "SBSA", "RNST", "TEST", "LTST", "DUNL", "ROSA", "PUSA", "PESA", "SPTS", "4587", "4588")) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNBT") {
    return(
      spp %>%
        dplyr::filter(species_code %in% c("BBWO", "ATTW")) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNOR") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Icterus") %>%
        pull(species_code)
    )
  }  
  if (annotation.code == "UNAM") {
    return(
      spp %>%
        dplyr::filter(species_code %in% c("AMRO", "SCTA", "WETA", "SUTA", "FCTA", "HETA", "BHGR", "RBGR")) %>%
        pull(species_code)
    )
  }  
  if (annotation.code == "UNKI") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Tyrannus") %>%
        pull(species_code)
    )
  }  
  if (annotation.code == "UNYE") {
    return(
      spp %>%
        dplyr::filter(species_code %in% c("GRYE", "LEYE")) |>
        pull(species_code)
    )
  }
  if (annotation.code == "USET") {
    return(
      spp %>%
        dplyr::filter(species_code %in% c("SPTO", "EATO")) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNPT") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Lagopus") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UMUR") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Uria") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UMUR") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Junco") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNVI") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Vireo") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNOW") {
    return(
      spp %>%
        dplyr::filter(species_order == "Strigiformes") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNSW") {
    return(
      spp %>%
        dplyr::filter(species_family == "Hirundinidae") %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UAHA") {
    return(
      spp %>%
        dplyr::filter(scientific_genus %in% c("Astur", "Accipiter", "Tachyspiza", "Lophospiza")) %>%
        pull(species_code)
    )
  } 
  if (annotation.code == "UNGR") {
    return(
      spp %>%
        dplyr::filter(species_family == "Phasianidae",
                      str_detect(
                        species_common_name,
                        regex(
                          "grouse|prairie-chicken|ptarmigan|capercaillie",
                          ignore_case = TRUE
                        )
                      )) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNRV") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Corvus",
                      str_detect(
                        species_common_name,
                        regex(
                          "raven",
                          ignore_case = TRUE
                        )
                      )) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNCR") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Corvus",
                      str_detect(
                        species_common_name,
                        regex(
                          "crow",
                          ignore_case = TRUE
                        )
                      )) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UJAE") {
    return(
      spp %>%
        dplyr::filter(scientific_genus == "Stercorarius",
                      str_detect(
                        species_common_name,
                        regex(
                          "jaeger",
                          ignore_case = TRUE
                        )
                      )) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNTW") {
    return(
      spp %>%
        dplyr::filter(scientific_genus %in% c("Pipilo", "Melozone")) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNBL") {
    return(
      spp %>%
        dplyr::filter(scientific_genus %in% c("Agelaius", "Quiscalus", "Euphagus", "Molothrus", "Xanthocephalus")) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UTEA") {
    return(
      spp %>%
        dplyr::filter(species_family == "Anatidae",
                      str_detect(
                        species_common_name,
                        regex(
                          "teal",
                          ignore_case = TRUE
                        )
                      )) %>%
        pull(species_code)
    )
  }
  if (annotation.code == "UNSS") {
    return(
      spp %>%
        dplyr::filter(species_code %in% c("NESP", "SALS")) %>%
        pull(species_code)
    )
  }

  #Use the more specific waterfowl and larid meanings where the unknown row
  #only supplies a shared family.
  if (annotation.code == "UNGO") {
    return(
      spp |>
        dplyr::filter(
          species_family == "Anatidae",
          str_detect(
            species_common_name,
            regex(
              "goose|brant",
              ignore_case = TRUE
            )
          )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "USWN") {
    return(
      spp |>
        dplyr::filter(scientific_genus == "Cygnus") |>
        pull(species_code)
    )
  }
  if (annotation.code == "UDAB") {
    return(
      spp |>
        dplyr::filter(
          scientific_genus %in%
            c(
              "Aix",
              "Anas",
              "Mareca",
              "Spatula"
            )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UDIV") {
    return(
      spp |>
        dplyr::filter(
          species_family == "Anatidae",
          !scientific_genus %in% c("Aix", "Anas", "Mareca", "Spatula", "Dendrocygna", "Cairina"),
          !str_detect(
            species_common_name,
            regex(
              "goose|brant|swan",
              ignore_case = TRUE
            )
          )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNME") {
    return(
      spp |>
        dplyr::filter(str_detect(
          species_common_name,
          regex(
            "merganser",
            ignore_case = TRUE
          )
        )) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNGU") {
    return(
      spp |>
        dplyr::filter(
          species_family == "Laridae",
          str_detect(
            species_common_name,
            regex(
              "gull|kittiwake",
              ignore_case = TRUE
            )
          )
        ) |>
        pull(species_code)
    )
  }
  if (annotation.code == "UNTE") {
    return(
      spp |>
        dplyr::filter(
          species_family == "Laridae",
          str_detect(
            species_common_name,
            regex(
              "\\btern\\b",
              ignore_case = TRUE
            )
          )
        ) |>
        pull(species_code)
    )
  }

  scientific.matches <- scientific.candidates(
    unknown.scientific.name,
    spp
  )
  if (length(scientific.matches)) {
    return(scientific.matches)
  }

  if (
    !is.na(unknown.family) &&
      unknown.family != "" &&
      unknown.family != "None"
  ) {
    return(
      spp |>
        dplyr::filter(species_family == unknown.family) |>
        pull(species_code)
    )
  }

  if (!is.na(unknown.order) && unknown.order != "") {
    return(
      spp |>
        dplyr::filter(species_order == unknown.order) |>
        pull(species_code)
    )
  }

  character()
}

#FORMAT WILDTRAX SPECIES LIST########

#1. Download the current WildTrax species object----
spp.wt <- wildrtrax::wt_get_species() |>
  mutate(across(where(is.character), str_squish))

#2. Format unknown bird categories from the WildTrax object----
unknown.categories <- spp.wt |>
  dplyr::filter(
    species_class == "AVES",
    str_detect(
      str_to_lower(species_common_name),
      "unknown|unidentified"
    )
  ) |>
  transmute(
    annotation_species_code = species_code,
    normalized_unknown_code = ifelse(
      species_code == "UNQK",
      "UNDU",
      ifelse(
        species_code == "UNAC",
        "UAHA",
        species_code
      )
    ),
    unknown_common_name = species_common_name,
    unknown_order = species_order,
    unknown_family = species_family,
    unknown_scientific_name = species_scientific_name,
    use_for_reclassification = !species_code %in%
      c(
        "UNBI",
        "178",
        "183"
      ),
    exclusion_reason = case_when(
      species_code == "UNBI" ~ "unbounded_unknown_bird",
      species_code %in% c("178", "183") ~ "legacy_numeric_code",
      TRUE ~ NA_character_
    )
  )

#UNKN is not included because WildTrax classifies it as a signal rather than
#an AVES species category.

#3. Format identified bird species from the same WildTrax object----
spp.known <- spp.wt |>
  dplyr::filter(
    species_class == "AVES",
    !species_code %in% unknown.categories$annotation_species_code,
    !is.na(species_scientific_name),
    species_scientific_name != ""
  ) |>
  mutate(
    scientific_name_key = str_to_upper(species_scientific_name),
    scientific_genus = word(species_scientific_name, 1)
  ) |>
  distinct(species_code, .keep_all = TRUE)

#4. Derive candidate species for every usable unknown category----
candidate.lookup <- unknown.categories |>
  dplyr::filter(use_for_reclassification) |>
  mutate(
    candidate_species_code = pmap(
      list(
        annotation_species_code,
        normalized_unknown_code,
        unknown_order,
        unknown_family,
        unknown_scientific_name
      ),
      candidate.species,
      spp = spp.known
    )
  ) |>
  dplyr::select(
    annotation_species_code,
    normalized_unknown_code,
    candidate_species_code
  ) |>
  unnest(candidate_species_code) |>
  left_join(
    spp.known |>
      dplyr::select(
        candidate_species_code = species_code,
        candidate_common_name = species_common_name,
        candidate_scientific_name = species_scientific_name,
        candidate_order = species_order,
        candidate_family = species_family
      ),
    by = "candidate_species_code"
  ) |>
  distinct()

unknown.categories <- unknown.categories |>
  left_join(
    candidate.lookup |>
      count(annotation_species_code, name = "number_candidate_species"),
    by = "annotation_species_code"
  ) |>
  mutate(
    number_candidate_species = replace_na(number_candidate_species, 0L),
    use_for_reclassification = use_for_reclassification &
      number_candidate_species > 0,
    exclusion_reason = case_when(
      number_candidate_species == 0 & is.na(exclusion_reason) ~
        "no_candidate_mapping",
      TRUE ~ exclusion_reason
    )
  )

#BUILD CLASSIFICATION DECISIONS######

#1. Set input and output database paths----
input.database <- if (test.run) {
  file.path(root, "WildTrax", "test", "01_wildtrax_raw_test.duckdb")
} else {
  file.path(
    root,
    "WildTrax",
    v.wt,
    paste0("01_wildtrax_raw_", v.wt, ".duckdb")
  )
}

output.database <- if (test.run) {
  file.path(
    root,
    "WildTrax",
    "test",
    "02b_unknown_classification_test.duckdb"
  )
} else {
  file.path(
    root,
    "WildTrax",
    v.wt,
    paste0("02b_unknown_classification_", v.wt, ".duckdb")
  )
}

dir.create(dirname(output.database), recursive = TRUE, showWarnings = FALSE)

#2. Save the lookup used for this run----
db.conn <- dbConnect(duckdb(), dbdir = output.database)
dbExecute(
  db.conn,
  paste0("ATTACH ", sql.path(input.database), " AS raw (READ_ONLY)")
)

#Remove tables written by earlier versions of this script. These objects are
#now either combined below or used only as temporary processing tables.
obsolete.tables <- c(
  "unknown_categories",
  "candidate_species_lookup",
  "unknown_classifier_scores",
  "unknown_classification_decisions",
  "unknown_reclassification_summary"
)
walk(
  obsolete.tables,
  ~ dbExecute(db.conn, paste("DROP TABLE IF EXISTS", .x))
)

unknown.species.lookup <- unknown.categories |>
  left_join(
    candidate.lookup,
    by = c(
      "annotation_species_code",
      "normalized_unknown_code"
    )
  )

dbWriteTable(
  db.conn,
  name = "unknown_species_lookup",
  value = unknown.species.lookup,
  overwrite = TRUE
)

#The two separate versions are retained temporarily because they make the SQL
#joins below simpler. They disappear when the database connection is closed.
dbWriteTable(
  db.conn,
  name = "unknown_categories_temp",
  value = unknown.categories,
  overwrite = TRUE,
  temporary = TRUE
)
dbWriteTable(
  db.conn,
  name = "candidate_species_lookup_temp",
  value = candidate.lookup,
  overwrite = TRUE,
  temporary = TRUE
)

#3. Select the unknown annotations----
#As in script 03, duplicate instances of a tag are removed before it is used.
dbExecute(
  db.conn,
  "
  CREATE OR REPLACE TEMP TABLE unknown_annotations AS
  SELECT
    a.project_id,
    a.task_id,
    a.recording_id,
    a.tag_id,
    a.species_code AS annotation_species_code,
    c.normalized_unknown_code,
    a.species_common_name AS annotation_common_name,
    c.use_for_reclassification,
    c.exclusion_reason,
    a.detection_time AS tag_start,
    COALESCE(a.tag_duration, 0) AS tag_duration,
    a.detection_time + COALESCE(a.tag_duration, 0) AS tag_end,
    a.abundance
  FROM raw.aru a
  INNER JOIN unknown_categories_temp c
    ON a.species_code = c.annotation_species_code
  QUALIFY ROW_NUMBER() OVER (
    PARTITION BY a.project_id, a.recording_id, a.tag_id
    ORDER BY a.detection_time
  ) = 1
  "
)

#4. Find and score eligible classifier species within the tag interval----
#One row is retained per candidate species and human annotation. The highest
#confidence is used if the classifier detected that species more than once.
dbExecute(
  db.conn,
  paste0(
    "
    CREATE OR REPLACE TEMP TABLE unknown_classifier_scores AS
    SELECT
      u.project_id,
      u.task_id,
      u.recording_id,
      u.tag_id,
      u.annotation_species_code,
      u.normalized_unknown_code,
      l.candidate_species_code,
      l.candidate_common_name,
      MAX(ai.confidence) AS classifier_score,
      COUNT(*) AS number_classifier_detections
    FROM unknown_annotations u
    INNER JOIN raw.aru_ai ai
      ON u.recording_id = ai.recording_id
      AND ai.detection_time >= u.tag_start - ",
    time.buffer.seconds,
    "
      AND ai.detection_time <= u.tag_end + ",
    time.buffer.seconds,
    "
    INNER JOIN candidate_species_lookup_temp l
      ON u.annotation_species_code = l.annotation_species_code
      AND ai.species_code = l.candidate_species_code
    WHERE ai.version = ",
    sql.string(classifier.version),
    "
      AND ai.is_species_allowed_in_project = TRUE
    GROUP BY ALL
    "
  )
)

#5. Rank the candidate species and build one decision per unknown tag----
#When only one candidate species is present, the runner-up score is zero.
#A tied top score consequently has a score difference of zero and is rejected.
dbExecute(
  db.conn,
  paste0(
    "
    CREATE OR REPLACE TABLE unknown_classification_results AS
    WITH ranked AS (
      SELECT
        *,
        ROW_NUMBER() OVER (
          PARTITION BY project_id, recording_id, tag_id
          ORDER BY classifier_score DESC, candidate_species_code
        ) AS score_rank
      FROM unknown_classifier_scores
    ),
    score_summary AS (
      SELECT
        project_id,
        recording_id,
        tag_id,
        MAX(CASE WHEN score_rank = 1 THEN candidate_species_code END)
          AS proposed_species_code,
        MAX(CASE WHEN score_rank = 1 THEN candidate_common_name END)
          AS proposed_common_name,
        MAX(CASE WHEN score_rank = 1 THEN classifier_score END)
          AS top_score,
        COALESCE(
          MAX(CASE WHEN score_rank = 2 THEN classifier_score END),
          0
        ) AS runner_up_score,
        COUNT(*) AS number_candidate_species
      FROM ranked
      GROUP BY project_id, recording_id, tag_id
    ),
    classifier_recordings AS (
      SELECT DISTINCT recording_id
      FROM raw.aru_ai
      WHERE version = ",
    sql.string(classifier.version),
    "
        AND is_species_allowed_in_project = TRUE
    )
    SELECT
      u.*,
      ",
    sql.string(classifier.version),
    " AS classifier_version,
      ",
    time.buffer.seconds,
    " AS time_buffer_seconds,
      ",
    minimum.score,
    " AS minimum_score,
      ",
    minimum.score.difference,
    " AS minimum_score_difference,
      s.proposed_species_code,
      s.proposed_common_name,
      s.top_score,
      s.runner_up_score,
      s.top_score - s.runner_up_score AS score_difference,
      COALESCE(s.number_candidate_species, 0) AS number_candidate_species,
      CASE
        WHEN u.use_for_reclassification = FALSE THEN
          COALESCE(u.exclusion_reason, 'not_eligible')
        WHEN r.recording_id IS NULL THEN 'no_classifier_report'
        WHEN s.proposed_species_code IS NULL THEN 'no_eligible_candidate'
        WHEN s.top_score <= ",
    minimum.score,
    " THEN 'score_below_threshold'
        WHEN s.top_score - s.runner_up_score <= ",
    minimum.score.difference,
    " THEN 'score_difference_below_threshold'
        ELSE 'accepted'
      END AS decision,
      CASE
        WHEN u.use_for_reclassification = TRUE
          AND s.top_score > ",
    minimum.score,
    "
          AND s.top_score - s.runner_up_score > ",
    minimum.score.difference,
    "
        THEN s.proposed_species_code
        ELSE NULL
      END AS reclassified_species_code
    FROM unknown_annotations u
    LEFT JOIN score_summary s
      ON u.project_id = s.project_id
      AND u.recording_id = s.recording_id
      AND u.tag_id = s.tag_id
    LEFT JOIN classifier_recordings r
      ON u.recording_id = r.recording_id
    "
  )
)

#6. Summarize unknowns before and after classification----
#Ignored and unmappable unknown types are included. Their number remaining is
#the same as their number before classification because they have no accepted
#reclassification decisions.
dbExecute(
  db.conn,
  "
  CREATE OR REPLACE TABLE unknown_classification_summary AS
  SELECT
    annotation_species_code,
    normalized_unknown_code,
    annotation_common_name,
    use_for_reclassification,
    exclusion_reason,
    COUNT(*) AS number_unknowns_before,
    COUNT(*) FILTER (WHERE decision = 'accepted') AS number_reclassified,
    COUNT(*) - COUNT(*) FILTER (WHERE decision = 'accepted')
      AS number_unknowns_after
  FROM unknown_classification_results
  GROUP BY ALL
  ORDER BY annotation_species_code
  "
)

#7. Print a short summary----
decision.summary <- dbGetQuery(
  db.conn,
  "
  SELECT decision, COUNT(*) AS number_annotations
  FROM unknown_classification_results
  GROUP BY decision
  ORDER BY number_annotations DESC
  "
)
print(decision.summary)

classification.results <- dbGetQuery(
  db.conn,
  "SELECT * FROM unknown_classification_summary"
)
print(classification.results)

#8. Format accepted reclassifications for review----
reclassified <- dbGetQuery(
  db.conn,
  "
  SELECT *
  FROM unknown_classification_results
  WHERE decision = 'accepted'
  ORDER BY annotation_species_code, reclassified_species_code
  "
)

reclassification.summary <- dbGetQuery(
  db.conn,
  "
  SELECT
    annotation_species_code,
    reclassified_species_code,
    proposed_common_name,
    COUNT(*) AS number_reclassified
  FROM unknown_classification_results
  WHERE decision = 'accepted'
  GROUP BY ALL
  ORDER BY annotation_species_code, number_reclassified DESC
  "
)
print(reclassification.summary)

message("Unknown classification decisions written to: ", output.database)
message("Test run: ", test.run)

dbExecute(db.conn, "DETACH raw")
dbDisconnect(db.conn, shutdown = TRUE)
