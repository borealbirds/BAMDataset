#!/bin/bash
# Submit from the BAMDataset repository with:
#   mkdir -p logs
#   sbatch 05.submit_covariates_slurm.sh
# Replace every value marked CHANGE ME before submission.

#SBATCH --account=def-group          # CHANGE ME: DRAC allocation
#SBATCH --job-name=bam_covariates
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100%20             # CHANGE ME: 1-number_of_manifest_rows
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

# CHANGE ME: load the R module used when BAMCovariates was installed.
module load StdEnv/2023 r/4.4.0

# CHANGE ME: paths visible from every compute node.
BAM_DATASET_REPOSITORY="/project/def-group/BAMDataset"
MANIFEST="/scratch/username/covariate_extraction/all_covariates/job_manifest.csv"

# Each task calls the same annotated controller used locally. The controller
# passes this Slurm array ID to BAMCovariates::run_job(). Valid completed jobs
# are skipped, so the same array may be resubmitted after an interruption.
Rscript --vanilla \
  "${BAM_DATASET_REPOSITORY}/05.CovariateExtraction.R" \
  drac_worker \
  "${MANIFEST}"
