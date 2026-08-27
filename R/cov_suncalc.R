# Covariate generating function for covariates based on sun position
#
# cov_name: character representing the desired name of the covariate
# type: which variable to calculate; either "altitude" (solar altitude angle), "azimuth" (solar azimuth angle), "d_altitude" (derivative of altitude), "d_azimuth" (derivative of azimuth)
# dt_seconds: numeric > 0; if type %in% c("d_altitude", "d_azimuth") - i.e., we are calculating a derivative, this argument determines the size of the finite difference approximation
#
# Returns a covariate generating function
cov_suncalc = function(cov_name,
                              type = c("altitude", "azimuth", "d_altitude", "d_azimuth"),
                              dt_seconds = 1) {
  
  library(suncalc)
  library(terra)
  library(lutz)
  library(sf)
  
  type = match.arg(type)
  
  deriv = str_detect(type, "^d_")
  var_keep = str_remove_all(type, "^d_")
  
  gfun = function(x, y, t, crs_in = add_EPSG(4326)) {
    
    # convert to degrees
    xy_proj = st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"), crs = crs_in) %>%
      st_transform("EPSG:4326")
    xy_proj_crds = st_coordinates(xy_proj)
    
    # get time zone associated with each time (we assume the times are in local time but need a time zone)
    tz_vals = tz_lookup(xy_proj, warn = FALSE)
    dates_forced = force_tzs(t, tz_vals)
    
    solar_data = data.frame(date = dates_forced,
                            lat = xy_proj_crds[, 2],
                            lon = xy_proj_crds[, 1])
    
    if (deriv) {
      solar_data_before = solar_data %>%
        mutate(date = date - dt_seconds)
      solar_data_after = solar_data %>%
        mutate(date = date + dt_seconds)
      
      spos_before = getSunlightPosition(data = solar_data_before, keep = var_keep)
      spos_after = getSunlightPosition(data = solar_data_after, keep = var_keep)
      
      return((spos_after[, 4] - spos_before[, 4]) / (dt_seconds * 2))
    }
    getSunlightPosition(data = solar_data, keep = var_keep)[, 4]
    
  }
  
  list(name = cov_name,
       domain = make_infinite_polygon(),
       get = gfun)
  
}