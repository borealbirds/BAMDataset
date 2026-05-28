# Generates the intersection of an aribtrary number of polygon objects, supplied as a list
#
# pl: SpatVector object or list of SpatVector objects
# crs_out: desired CRS of output object (by default, selects the first polygon)
# do_sf: if TRUE, converts everything to sf before proceeding
#
# Returns a SpatVector object representing the intersection of all polygons
intersect_many = function(pl,
                          crs_out,
                          do_sf = FALSE) {
  
  # function that only projects if the CRS actually needs to change - I don't believe terra does this by default - and hopefully speeds things up a bit
  projfun = function(x) {
    if (is_infinite_polygon(x)) return(make_infinite_polygon(crs_out))
    if (crs(crs_out) != crs(x)) x = terra::project(x, crs_out)
    x
  }
  
  # if it's not a list and just one object, just return the object
  if (is(pl, "SpatVector")) pl = list(pl)
  if (missing(crs_out)) crs_out = crs(pl[[1]])
  if (length(pl) == 1) return(projfun(pl[[1]]))
  
  pl_proj = lapply(pl, projfun)
  if (do_sf) {
    pl_proj = lapply(pl_proj, st_as_sf)
    sf_int = st_intersection(do.call(st_as_sfc, lds))
    return(vect(sf_int))
  }
  
  out = pl_proj[[1]] * pl_proj[[2]]
  
  for (i in seq_len(length(pl_proj) - 2))  out = out * pl_proj[[i + 2]]
  
  out
  
}

# helper function - should be enough for most polygons, hopefully?
is_infinite_polygon = function(x) {
  if (nrow(crds(x)) > 5) return(FALSE)
  all(crds(make_infinite_polygon(crs(x))) == crds(x))
}
