# Covariate generating function for covariates that are temporally and spatially variable
#
# cov_name: character representing the desired name of the covariate
# raster_dir: folder where all the raster files are located. the files inside will be named according to dates (and not have any other content in the names)
# by: character representing POSIXct format expected in file names (e.g., the default option of "%Y" implies that the file names will all include years).
# method: character - either "nearest" (matches the data to the nearest available date) or "exact" (dates must match exactly or else an error is thrown)
# buffer: numeric > 0; buffers points by a certain amount and extracts a function (see ...) of those values
# ...: additional arguments to exact_extract, only relevant if buffer > 0
#
# Returns a covariate generating function
cov_dynamic_raster = function(cov_name,
                              raster_dir,
                              by = "%Y",
                              method = c("exact", "nearest"), 
                              buffer = 0,
                              ...) {
  
  library(terra)
  library(sf)
  library(exactextractr)
  
  method = match.arg(method)
  
  raster_files = list.files(raster_dir, full.names = TRUE, pattern = "\\.tif")
  raster_files_dates_only = list.files(raster_dir, pattern = "\\.tif") %>% str_split_i(".tif", 1)
  raster_files_dates_pos = as.POSIXct(raster_files_dates_only, format = by)
  
  all_rasters = lapply(raster_files, rast)
  all_crs_values = sapply(all_rasters, crs)
  all_file_polygons = lapply(all_rasters, function(r) vect(ext(r), crs = crs(r)))
  
  gfun = function(x, y, t, crs_in = all_crs_values[1]) {
    vdf = data.frame(x = x, y = y, t = t) %>%
      mutate(tf = as.character(format(t, by)),
             ord = 1:n(),
             file_ind = match(tf, raster_files_dates_only)) %>%
      arrange(file_ind)
    
    if (any(is.na(vdf$file_ind))) {
      if (method == "exact") stop("Not all input values have matching covariate rasters. Make sure 'raster_dir' has all necessary files and that 'by' is properly specified for the naming conventions of each raster file.")
      
      vdf_na = vdf %>% 
        dplyr::filter(is.na(file_ind)) %>%
        mutate(file_ind = sapply(as.POSIXct(tf, format = by), function(dd) which.min(abs(dd - raster_files_dates_pos))))
      
      vdf[is.na(vdf$file_ind), ] = vdf_na
      vdf = vdf %>% arrange(file_ind, ord)
    }
    
    vdf = vdf %>% mutate(crs_out = all_crs_values[file_ind]) # do this at the end so any NA's are fixed
    
    if (buffer > 0) {
      # make buffered polygons
      v_sf = st_as_sf(vdf, coords = c("x", "y"), crs = crs_in)
      v_buf = st_buffer(v_sf, buffer)
    }
    
    vdf$cov = unlist(lapply(unique(vdf$file_ind), function(ind) {
      
      this_file = vdf[vdf$file_ind == ind, ]
      if (buffer > 0) return(exact_extract(all_rasters[[ind]], st_transform(v_buf[vdf$file_ind == ind, ], this_file$crs_out[1]), ...))
      
      this_file_vect = vect(this_file, geom = c("x", "y"), crs = crs_in) %>% terra::project(this_file$crs_out[1])
      terra::extract(all_rasters[[ind]], this_file_vect, ID = FALSE)[, 1]
      
    }))
    
    vdf = vdf %>% 
      arrange(ord) %>% # bring back to original order
      mutate(cov = ifelse(is.na(cov), 0, cov))
    vdf$cov # we only return the covariate values, nothing else
    
  }
  
  list(name = cov_name,
       domain = intersect_many(all_file_polygons),
       get = gfun)
  
}