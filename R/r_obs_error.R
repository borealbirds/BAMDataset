# Helper function for path simulation; adds measurement / observation error to movement paths
#
# dt_in: data.table with columns representing information about 
# sigma_o: observation error variance (if equal to 0, no spatial error included)
# obs_keep: number between 0 and 1 representing the percentage of observations to be retained
# id_keep: number between 0 and 1 representing the percentage of observations for which individual ID has been retained
# bound: geom_sf polygon object representing areas where animals cannot travel. When we extract values, should give a non-NA value for eligible locations.
# x_scale: the length of one spatial unit, in the units associated with the movement parameters
# range_centre: matrix with 2 columns and the same number of rows as prev_state, indicating (unscaled) range centres
#
# Returns a data.table with new columns to reflect measurement error. Note that "true" locations (as well as locations that were "removed") are kept in case they're needed.
r_obs_error = function(dt_in,
                       sigma_o = 0,
                       obs_keep = 1,
                       id_keep = 1,
                       bound = NULL,
                       x_scale = 1,
                       range_centre = dt_in[, .("x", "y")] * 0) {
  
  dt_in_witherrors = dt_in[, c("x_obs", "y_obs") := .(x + rnorm(.N, sd = sigma_o), y + rnorm(.N, sd = sigma_o))][ # Jitter based on observation error
    , obs_kept := rbinom(.N, 1, obs_keep)][ # Remove some observations entirely
      , id_kept := ifelse(rbinom(.N, 1, id_keep), animal_id, NA)] # remove some ID's entirely
  
  if (is.null(bound)) return(dt_in_witherrors)
  
  # Fix any jitters that jittered into the boundary
  locs_sf = st_as_sf(dt_in_witherrors[, .(x_obs, y_obs)] * x_scale + range_centre, coords = c("x_obs", "y_obs"), crs = st_crs(bound))
  na_ext_vals = !(1:nrow(dt_in) %in% st_intersects(bound, locs_sf)[[1]])
  
  if (!any(na_ext_vals)) return(dt_in_witherrors)
  
  n_resim = sum(na_ext_vals)
  # Get nearest points within the boundary for each location that's bad
  nearest_inbounds_pts = st_nearest_points(locs_sf[na_ext_vals, ], bound)
  nearest_pts_coords = st_coordinates(nearest_inbounds_pts)[seq(from = 2, to = n_resim * 2, by = 2), 1:2, drop = FALSE]
  # shift back to "model" units and coordinates
  nearest_pts_shifted = (nearest_pts_coords - range_centre[na_ext_vals, ]) / x_scale
  # Add in new values
  dt_in_witherrors$x_obs[na_ext_vals] = nearest_pts_shifted[, 1]
  dt_in_witherrors$y_obs[na_ext_vals] = nearest_pts_shifted[, 2]
  
  dt_in_witherrors
  
}