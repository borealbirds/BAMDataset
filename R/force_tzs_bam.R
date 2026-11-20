# Temporarily force local clock times to their coordinate-based time zones
# using base R.
#
# This helper exists because timechange 0.4.0 cannot read the updated TZif
# files for several Canadian time zones in the tzdata 2026e files bundled with
# R 4.6.2. This includes America/Edmonton, America/Vancouver,
# America/Inuvik, and America/Winnipeg. Base R reads these files correctly and
# retains their historical time-zone rules.
#
# Check this compatibility issue when timechange is updated. Remove this helper
# and return to lubridate::force_tzs() in cov_timeofday.R and cov_suncalc.R once
# timechange:::C_valid_tz("America/Edmonton") returns TRUE with current tzdata.

force_tzs_bam = function(t, tz) {
  time_text = format(
    t,
    format = "%Y-%m-%d %H:%M:%OS",
    tz = "UTC",
    usetz = FALSE
  )
  out = as.POSIXct(rep(NA_character_, length(t)), tz = "UTC")

  for (this_tz in unique(tz[!is.na(tz)])) {
    use = tz == this_tz & !is.na(t)
    out[use] = as.POSIXct(
      time_text[use],
      format = "%Y-%m-%d %H:%M:%OS",
      tz = this_tz
    )
  }

  out
}

# suncalc::getSunlightTimes() also calls timechange internally using the
# computer's time zone. Temporarily use UTC for that internal calculation and
# restore the user's setting immediately afterward. Remove this wrapper along
# with force_tzs_bam once timechange compatibility is restored.
get_sunlight_times_bam = function(...) {
  old_tz = Sys.getenv("TZ", unset = NA_character_)
  Sys.setenv(TZ = "UTC")
  on.exit({
    if (is.na(old_tz)) {
      Sys.unsetenv("TZ")
    } else {
      Sys.setenv(TZ = old_tz)
    }
  })

  suncalc::getSunlightTimes(...)
}
