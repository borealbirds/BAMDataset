# Covariate generating function for covariates based on sun position
#
# cov_name: character representing the desired name of the covariate
# type: which variable to calculate; either "sunrise" for time since sunrise or "sunset" for time since sunset
# units: for "difftime", what units to return for time since
#
# Returns a covariate generating function
cov_timeofday = function(cov_name,
                         type = c("sunrise", "sunset"),
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