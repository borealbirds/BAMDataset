# BAMDataset

Code for downloading and updating data that BAM has access to from WildTrax

## Covariate extraction

`05.CovariateExtraction.R` is the single user-facing controller for extracting
survey covariates with the `BAMCovariates` package. Its settings block supports
four modes:

- `local`: build restartable survey-chunk jobs, run them locally, combine the
  results, and optionally update the existing BAMDataset DuckDB;
- `drac_prepare`: create the manifest and survey chunks on DRAC without running
  extraction;
- `drac_worker`: run one manifest row inside a Slurm array; and
- `combine`: validate completed chunks, combine them, and optionally update the
  DuckDB after the Slurm array finishes.

Start with `test.run <- TRUE` and `write.covariate.table <- FALSE`. Inspect the
resulting `covariates` object before allowing the production database update.
Existing covariate columns are retained by default. Set
`replace.existing.outputs <- TRUE` only when deliberately recalculating the
outputs selected by `stage.outputs`.

For a DRAC run:

1. Stage the DuckDB and required native rasters on storage visible to compute
   nodes.
2. Set `BAM_DATASET_ROOT`, `BAM_SPATIAL_ROOT`,
   `BAM_COVARIATE_JOB_ROOT`, and `BAM_DATASET_REPOSITORY`.
3. Run `05.CovariateExtraction.R` with `workflow.mode <- "drac_prepare"`.
4. Edit the allocation, manifest path, and array range in
   `05.submit_covariates_slurm.sh`, then submit it with `sbatch`.
5. After all workers succeed, run the controller in `combine` mode and inspect
   `covariates` before enabling the DuckDB write.

When DRAC cannot store every raster at once, assign a unique `stage.name` and
list the available recipe `output_name` values in `stage.outputs`. Each stage
uses its own checkpoint directory. Job and combined RDS files are restart
checkpoints; the DuckDB `covariate` table remains the final dataset product.
