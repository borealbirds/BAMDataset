# Covariate generating function for covariates based on sun position
#
# cov_name: character representing the desired name of the covariate
# type: which variable to calculate; either "sunrise" for time since sunrise or "sunset" for time since sunset
# units: for "difftime", what units to return for time since
#
# Returns a covariate generating function
cov_timeofday = function(cov_name,
                         type = c("sunrise", "solarNoon", "nadir", "sunset", "sunriseEnd", "sunsetStart", "dawn", "dusk", "nauticalDawn", "nauticalDusk", "nightEnd", "night", "goldenHourEnd", "goldenHour"),
                         units = c("hours", "secs", "mins", "days", "weeks")) {
  
  library(suncalc)
  library(terra)
  library(lutz)
  library(sf)
  
  type = match.arg(type)
  units = match.arg(units)
  
  gfun = function(x, y, t, crs_in = add_EPSG(4326)) {
    
    # convert to degrees
    xy_proj = st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = crs_in) %>%
      st_transform("EPSG:4326")
    xy_proj_crds = st_coordinates(xy_proj)
    
    # get time zone associated with each time (we assume the times are in local time but need a time zone)
    tz_vals = tz_lookup(xy_proj, warn = FALSE)
    dates_forced = force_tzs(t, tz_vals)
    
    solar_data = data.frame(date = as.Date(t),
                            lat = xy_proj_crds[, 2],
                            lon = xy_proj_crds[, 1])
    
    sun_times = getSunlightTimes(data = solar_data, keep = type)[, 4]
    as.numeric(difftime(dates_forced, sun_times, units = units))
    
  }
  
  list(name = cov_name,
       domain = make_infinite_polygon(),
       get = gfun)
  
}

# Get discrete time of day classification
#
# t_since_dawn: numeric; # of hours after beginning of morning civil twilight
# t_since_sunrise: numeric; # of hours after sunrise
# t_since_goldenend: numeric; # of hours after the end of morning golden hour
# t_since_golden: numeric; # of hours after beginning of evening golden hour
# t_since_dusk: numeric; # of hours after beginning of evening nautical twilight
# t_since_nadir: numeric; # of hours after the darkest moment of the night / day
#
# Returns a character vector of values
get_tod = function(t_since_dawn, 
                   t_since_sunrise, 
                   t_since_goldenend, 
                   t_since_golden, 
                   t_since_dusk, 
                   t_since_nadir) {
  
  case_when(
    # 1. sunrise window: within golden hour
    !is.na(t_since_dawn) & t_since_dawn >= 0 & t_since_goldenend <= 0 ~ "sunrise",
    # 2. sunrise window: when there is no dawn anything from nadir to end of goldenhour
    is.na(t_since_dawn) & !is.na(t_since_sunrise) & t_since_nadir >= 0 & t_since_goldenend <= 0 ~ "sunrise",
    # 3. sunrise window: when there is no sunrise and then time is during golden hours after nadir
    is.na(t_since_sunrise) & t_since_nadir >= 0 & t_since_goldenend <= 0 ~ "sunrise",
    # 4. dawn: between goldenhour and +5h after sunrise
    !is.na(t_since_sunrise) & t_since_goldenend > 0 & t_since_sunrise <= 5 ~ "dawn",
    # 5. dawn: define dawn when there is no sunrise
    is.na(t_since_dawn) & is.na(t_since_sunrise) & t_since_nadir >= 0 & t_since_nadir <= 5 ~ "dawn",
    # 6. sunset window: between evening goldenhour and dusk when dusk is before midnight
    !is.na(t_since_dusk) & t_since_golden > t_since_dusk & t_since_golden >= 0 & t_since_dusk <= 0 ~ "sunset",
    # 7. sunset window: between evening goldenhour and dusk when dusk is after midnight
    !is.na(t_since_dusk) & t_since_golden < t_since_dusk & (t_since_golden >= 0 | t_since_dusk <= 0) ~ "sunset",
    # 8. sunset window: when there is no dusk anything from goldenhour to nadir; nadir is before midnight
    is.na(t_since_dusk) & t_since_golden > t_since_nadir & t_since_golden >= 0 & t_since_nadir < 0 ~ "sunset",
    # 9. sunset window: when there is no dusk anything from goldenhour to nadir; nadir is after midnight
    is.na(t_since_dusk) & t_since_golden < t_since_nadir & (t_since_golden >= 0 | t_since_nadir < 0) ~ "sunset",
    # 10. night: between dusk and dawn, regardless of whether dusk falls before or after midnight
    # if dusk comes after midnight (ie. dusk < dawn) # if dusk comes before midnight (ie. dusk > dawn)
    !is.na(t_since_dusk) & !is.na(t_since_dawn) & if_else(t_since_dawn < t_since_dusk, t_since_dusk >= 0 & t_since_dawn <= 0, t_since_dusk >= 0 | t_since_dawn <= 0) ~ "night",
    # 11. For all other cases it should be day
    TRUE ~ "day"
  )
  
}
