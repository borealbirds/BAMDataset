# Helper function for above
make_infinite_polygon = function(crs = add_EPSG(4326), return_type = c("terra", "sf"), nseq = 10) {
  
  return_type = match.arg(return_type)
  
  x_seq = seq(-180, 180, length.out = nseq)
  y_seq = seq(-90, 90, length.out = nseq)
  xy = cbind(rep(x_seq, length(y_seq)), rep(y_seq, each = length(x_seq)))
  v_bounds = vect(xy, crs = add_EPSG(4326)) # this is always in lat-lon
  v_bounds_proj = terra::project(v_bounds, crs)
  v_ext = as.polygons(ext(v_bounds_proj), crs = crs)
  if (return_type == "sf") v_ext = st_as_sf(v_ext)
  v_ext
  
}