# Helper function for path simulation; adds spacing information to data.table, incorporating boundaries if necessary
#
# dt_in: data.table with columns representing information about centroid location but not (yet) about spacing
# ssx_means: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): mean individual spacing distances in x direction for each row of dt_in
# ssy_means: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): mean individual spacing distances in y direction for each row of dt_in
# sfx_means: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): mean foraging spacing distances in x direction for each row of dt_in
# sfy_means: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): mean foraging spacing distances in y direction for each row of dt_in
# ss_sds: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): individual spacing standard deviations for each row of dt_in
# sy_sds: vector of length equal to # of rows in dt_in (or 1, in which case the same value will be repeated): foraging spacing standard deviations for each row of dt_in
# bound: geom_sf polygon object representing areas where animals cannot travel. When we extract values, should give a non-NA value for eligible locations.
# W: SpatRaster with information about habitat quality. Values should all be greater than or equal to 0.
# n_random: integer > 0; number of random draws to generate for each group. If habitat is not defined (i.e. NULL) this value will be set to 1.
# x_scale: the length of one spatial unit, in the units associated with the movement parameters
# range_centre: matrix with 2 columns and the same number of rows as prev_state, indicating (unscaled) range centres
# 
# Returns a data.table with added columns to reflect spacing
r_spacing = function(dt_in,
                     ssx_means = 0,
                     ssy_means = 0,
                     sfx_means = 0,
                     sfy_means = 0,
                     ss_sds = 0,
                     sf_sds = 0,
                     bound = NULL,
                     W = NULL,
                     n_random = 1,
                     x_scale = 1,
                     range_centre = dt_in[, .("x", "y")] * 0) {
  
  if (is.null(W) && n_random != 1) {
    warning("n_random should be set to 1 if W (habitat raster) is NULL. Correcting now...")
    n_random = 1
  }
  
  n_data = nrow(dt_in)
  
  # Figure out how many times we need to repeat each value to make it work
  rp_fac_ssx_means = (length(ssx_means) == 1) * n_random * (n_data - 1) + n_random
  rp_fac_ssy_means = (length(ssy_means) == 1) * n_random * (n_data - 1) + n_random
  rp_fac_sfx_means = (length(sfx_means) == 1) * n_random * (n_data - 1) + n_random
  rp_fac_ssy_means = (length(sfy_means) == 1) * n_random * (n_data - 1) + n_random
  rp_fac_ss_sds = (length(ss_sds) == 1) * n_random * (n_data - 1) + n_random
  rp_fac_sf_sds = (length(sf_sds) == 1) * n_random * (n_data - 1) + n_random
  
  # Repeat them
  ssx_means = rep(ssx_means, each = rp_fac_ssx_means)
  ssy_means = rep(ssy_means, each = rp_fac_ssy_means)
  sfx_means = rep(sfx_means, each = rp_fac_sfx_means)
  sfy_means = rep(sfy_means, each = rp_fac_ssy_means)
  ss_sds = rep(ss_sds, each = rp_fac_ss_sds)
  sf_sds = rep(sf_sds, each = rp_fac_sf_sds)
  
  dt_withspacing = dt_in[rep(1:.N, each = n_random)][
    , c("ssx", "ssy", "sfx", "sfy") := .(rnorm(.N, ssx_means, ss_sds),
                                         rnorm(.N, ssy_means, ss_sds),
                                         rnorm(.N, sfx_means, sf_sds),
                                         rnorm(.N, sfy_means, sf_sds))][
                                           , c("x", "y") := .(cpx + ssx + sfx,
                                                              cpy + ssy + sfy)]
  
  if (!is.null(bound)) {
    # Correct any values that land within the boundaries
    locs_sf = st_as_sf(dt_withspacing[, .(x, y)] * x_scale + range_centre, coords = c("x", "y"), crs = st_crs(bound))
    na_ext_vals = !(1:nrow(dt_in) %in% st_intersects(bound, locs_sf)[[1]])
    
    if (any(na_ext_vals)) {
      n_resim = sum(na_ext_vals)
      # Get nearest points within the boundary for each location that's bad
      nearest_inbounds_pts = st_nearest_points(locs_sf[na_ext_vals, ], bound)
      nearest_pts_coords = st_coordinates(nearest_inbounds_pts)[seq(from = 2, to = n_resim * 2, by = 2), 1:2, drop = FALSE]
      # shift back to "model" units and coordinates
      nearest_pts_shifted = (nearest_pts_coords - range_centre[na_ext_vals, ]) / x_scale
      # Add in new values
      tot_sp_ratios = (nearest_pts_shifted - dt_withspacing[na_ext_vals, .(cpx, cpy)]) / (dt_withspacing[na_ext_vals, .(ssx, ssy)] + dt_withspacing[na_ext_vals, .(sfx, sfy)])
      dt_withspacing$x[na_ext_vals] = nearest_pts_shifted[, 1]
      dt_withspacing$y[na_ext_vals] = nearest_pts_shifted[, 2]
      # Re-calculate spacing values
      dt_withspacing[na_ext_vals, c("ssx", "ssy", "sfx", "sfy")] = dt_withspacing[na_ext_vals, .(ssx, ssy, sfx, sfy)] * cbind(tot_sp_ratios, tot_sp_ratios)
    }
  }
  
  if (!is.null(W)) {
    # Do the habitat
    res_coords = dt_withspacing[, .(x, y)] * x_scale + range_centre
    dt_withspacing$habitat = terra::extract(W, res_coords, ID = FALSE)
    # Correct NA values to 0 if necessary.
    dt_withspacing$habitat[is.na(dt_withspacing$habitat)] = 0
    # If all habitat values are 0 set them to 1, for now. This may have unintended consequences but we'll get to that later.
    dt_withspacing = dt_withspacing[, sum_habitat := sum(habitat), by = path_id][
      , habitat := ifelse(sum_habitat == 0, 1, habitat)]
    # Sample randomly
    dt_withspacing = dt_withspacing[, .SD[sample(.N, 1, prob = habitat)], by = .(path_id, animal_id)]
  }
  
  dt_withspacing
  
}