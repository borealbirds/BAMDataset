# Bootstrapped log offsets for many surveys at once, using matrices instead of a long data.frame
#
# Replaces the "pars_df" / get_offset_la_mc / group_by / summarize block in 04.CalculateOffsets.R.
# Linear predictors are built as (surveys x parameters) %*% (parameters x bootstrap draws), offsets are
# computed with get_offset_la_cf, and the mean and SD across draws are taken per survey. Surveys are
# processed in chunks so memory stays bounded even for millions of surveys.
#
# dat: data.frame of surveys with columns timeofday, timeofyear, biomass, near_road, duration_min, distance_km
# r_pars: data.frame / matrix of random parameter draws (one row per bootstrap draw, named columns)
# include_road: logical; should the near_road effect on log_alpha be used for this species?
# tod_to_keep: character vector of time-of-day bins for which offsets are reported (others become NA)
# all_par_names: all parameter names that may appear; missing ones are filled with 0 (as check_names does)
# chunk_size: number of surveys to process at once
#
# Returns a data.frame with columns mean_log_offset and sd_log_offset, one row per row of dat
get_offset_boot_mat = function(dat,
                               r_pars,
                               include_road = FALSE,
                               tod_to_keep = NULL,
                               all_par_names = character(0),
                               chunk_size = 2e5) {

  tod_levels = c("nauticaldawn", "sunrise", "morning", "day", "sunset", "nauticaldusk", "night")

  r_pars = check_names(as.data.frame(r_pars), all_par_names, fill = 0)
  r_pars = check_names(r_pars, c(paste0("log_lambda_", tod_levels), "log_lambda_timeofyear",
                                 "log_alpha_", "log_alpha_near_road", "log_alpha_open_closed"), fill = 0)
  n_boot = nrow(r_pars)

  # parameters x draws
  B_lambda = t(as.matrix(r_pars[, c(paste0("log_lambda_", tod_levels), "log_lambda_timeofyear")]))
  B_alpha = t(as.matrix(r_pars[, c("log_alpha_", "log_alpha_near_road", "log_alpha_open_closed")]))

  # surveys x parameters
  tod = as.character(dat$timeofday)
  X_lambda = cbind(outer(tod, tod_levels, "==") * 1, dat$timeofyear)
  use_road = !is.na(include_road) && include_road
  X_alpha = cbind(1, if (use_road) dat$near_road else 0, dat$biomass)

  n = nrow(dat)
  mean_log_offset = sd_log_offset = rep(NA_real_, n)
  keep = tod %in% tod_to_keep

  for (start in seq(1, n, by = chunk_size)) {
    idx = start:min(n, start + chunk_size - 1)
    idx = idx[keep[idx]]
    if (length(idx) == 0) next

    log_lambda = X_lambda[idx, , drop = FALSE] %*% B_lambda
    log_alpha = X_alpha[idx, , drop = FALSE] %*% B_alpha

    # column-major: each block of length(idx) values is one bootstrap draw
    off = matrix(get_offset_la_cf(c(log_lambda), c(log_alpha),
                                  rep(dat$duration_min[idx], n_boot),
                                  rep(dat$distance_km[idx], n_boot)),
                 ncol = n_boot)

    m = rowMeans(off)
    mean_log_offset[idx] = m
    sd_log_offset[idx] = if (n_boot > 1) sqrt(rowSums((off - m) ^ 2) / (n_boot - 1)) else NA_real_
  }

  data.frame(mean_log_offset = mean_log_offset,
             sd_log_offset = sd_log_offset)

}
