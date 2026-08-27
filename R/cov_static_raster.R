# Covariate generating function for covariates that are temporally static and spatially variable
#
# cov_name: character representing the desired name of the covariate
# raster_path: character representing the file path (including file name and extension) of the raster containing information about the covariates
# fun_app: function to apply to results
# cov_scale: numeric > 0; scale factor to apply to results BEFORE fun_app is called
#
# Returns a covariate generating object with three elements: "name" (the name of the covariate), "domain" (SpatVector showing where the covariate is defined), and "get" (function for getting covariate values)
cov_static_raster = function(cov_name, 
                             raster_path,
                             fun_app = function(x) {x},
                             cov_scale = 1) {
  
  r = rast(raster_path)
  
  gfun = function(x, y, t, crs_in = crs(r)) {
    v = vect(cbind(x, y), crs = crs_in) %>% terra::project(crs(r))
    fun_app(cov_scale * terra::extract(r, v, ID = FALSE)[, 1])
  }
  
  list(name = cov_name,
       domain = as.polygons(r * 0, values = FALSE),
       get = gfun)
  
}
