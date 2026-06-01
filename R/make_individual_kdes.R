# Generate a polygon (or multipolygon) of MCP home ranges for a series of animal tracks (useful for pre-processing and viewing data)
#
# x: numeric vector of x coordinates (longitude or UTM)
# y: numeric vector of y coordinates (latitude or UTM)
# t: vector of dates and times as POSIXct
# id: character vector with individual identities. If not supplied, assumes there is only one animal in the dataset
# crs_in: integer representing EPSG code associated with CRS of input data. Defaults to 4326 (WGS84 lat/lon), so if your input data are not in degrees, please change.
# crs_desired: integer representing EPSG code associated with CRS at which analyses will be done. Defaults to crs_in.
# verbose: integer >= 0; higher values represent more print output during the course of function evaluation; currently only coded up to 1
# buffer: numeric >= 0; spatial buffer around polygons in m
# min_locs: integer > 0; 
# ...: additional arguments to hr_isopleths (e.g., levels)
# 
# Returns a set of SpatVector polygons including named ID and each representing a MCP home range
make_individual_kdes = function(x,
                                y,
                                t,
                                id = character(length(x)),
                                crs_in = 4326,
                                crs_desired = crs_in,
                                verbose = 0,
                                buffer = 0,
                                min_locs = 1,
                                ...) {
  
  require(amt)
  
  if (crs_in == 4326 & (min(x) < -180 | min(y) < -90 | max(x) > 180 | max(y) > 90)) warning("Some coordinate values (x or y) appear to be outside the expected bounds for latitudes and longitudes, yet input CRS is set to EPSG:4326 (the default choice). Please confirm that this is the correct value for this argument or supply a different value.")
  
  all_ids = unique(id)
  
  xy = cbind(x, y)
  if (crs_desired != crs_in) xy = terra::project(xy, add_EPSG(crs_in), add_EPSG(crs_desired))
  
  poly_out_list = lapply(all_ids, function(ind) {
    if (verbose > 0) message("Beginning ID ", ind)
    id_bool = id == ind # which values come from this ID?
    if (sum(id_bool) < min_locs) return(vect()) # return blank polygon here because weird errors begin to pop up with 4 or fewer locations
    xy_track = make_track(data.frame(x = xy[id_bool, 1], y = xy[id_bool, 2], t = t[id_bool]), .x = x, .y = y, .t = t)
    this_kde_poly = vect(hr_isopleths(hr_kde(xy_track), ...))
    crs(this_kde_poly) = add_EPSG(crs_desired)
    this_kde_poly
  })
  
  all_ids = all_ids[sapply(poly_out_list, nrow) > 0]
  poly_out = do.call(rbind, poly_out_list)
  poly_out$name = all_ids # so polygons keep ID names from the input data
  if (buffer > 0) poly_out = buffer(poly_out, buffer)
  poly_out
  
}
