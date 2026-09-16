# Helper function for path simulation; generates random iterates from the Ornstein-Uhlenbeck with Foraging (OUF) model given a previous set of state values
#
# prev_state: matrix or data.frame with 4 columns indicating the previous values of cx (centroid x), cy (centroid y), cvx (velocity in x direction), and cvy (velocity in y direction) for each component **in that order**
# H: state transition matrix (2x2)
# R: state covariance matrix (2x2)
# sigma_v: variability in velocity
# beta_v: temporal autocorrelation associated with velocity
# omega_v: intensity of range restriction
# bound: geom_sf polygon object representing areas where animals cannot travel. When we extract values, should give a non-NA value for eligible locations.
# W: SpatRaster with information about habitat quality. Values should all be greater than or equal to 0.
# n_random: integer > 0; number of random draws to generate for each group. If habitat is not defined (i.e. NULL) this value will be set to 1.
# do_repeat: integer > 0; includes the centroid more than once in the data.table (often useful for groups of more than one animal)
# t_label: value to include for the "time" column of the data.table
# x_scale: the length of one spatial unit, in the units associated with the movement parameters
# t_scale: the length of one timestep, in the units associated with the movement parameters, only relevant if H and R are not both supplied
# range_centre: matrix with 2 columns and the same number of rows as prev_state, indicating (unscaled) range centres
# ...: additional arguments to correct_direction_bound (direction adjustment algorithm)
#
# Produces a data.frame with new information
r_ouf = function(prev_state,
                 H,
                 R,
                 sigma_v,
                 beta_v,
                 omega_v,
                 bound = NULL,
                 W = NULL,
                 n_random = 1,
                 do_repeat = 1,
                 t_label = 0,
                 x_scale = 1,
                 t_scale = 1,
                 range_centre = prev_state[, 1:2, drop = FALSE] * 0,
                 ...) {
  
  require(data.table)
  
  if (missing(H) || missing(R)) {
    # make H and R from sigma_v, beta_v, omega_v
    if (missing(sigma_v) || missing(beta_v) || missing(omega_v)) stop("must either supply H and R, or sigma_v, beta_v, omega_v")
    
    omega2_v = omega_v ^ 2
    nu_v = sqrt(beta_v ^ 2 - omega2_v)
    
    if (nu_v == 0) {
      C0 = exp(-t_scale * beta_v) * (1 + t_scale * beta_v)
      C1 = -omega2_v * t_scale * exp(-t_scale * beta_v)
      C2 = omega2_v * exp(-t_scale * beta_v) * (1 - t_scale * beta_v)
    } else {
      C0 = exp(t_scale * (nu_v - beta_v)) * (1 + beta_v / nu_v) / 2 + exp(t_scale * (-nu_v - beta_v)) * (1 - beta_v / nu_v) / 2
      C1 = -omega2_v / nu_v / 2 * (exp(t_scale * (nu_v - beta_v)) - exp(t_scale * (-nu_v - beta_v)))
      C2 = -omega2_v / 2 * (exp(t_scale * (nu_v - beta_v)) * (1 - beta_v / nu_v) + exp(t_scale * (-nu_v - beta_v)) * (1 + beta_v / nu_v))
    }
    
    V_pp = sigma_v * (1 - C0 ^ 2 - C1 ^ 2 / omega2_v)
    V_vp = 2 * beta_v * sigma_v * C1 ^ 2 / omega2_v
    V_vv = sigma_v * (omega2_v - C1 ^ 2 - C2 ^ 2 / omega2_v)
    H_pp = C0
    H_pv = C1
    H_vp = -C1 / omega2_v
    H_vv = -C2 / omega2_v
    
    H = matrix(nrow = 2, ncol = 2, data = c(H_pp, H_pv, H_vp, H_vv))
    R = matrix(nrow = 2, ncol = 2, data = c(V_pp, V_vp, V_vp, V_vv))
    
  }
  
  if (is.null(W) && n_random != 1) {
    warning("n_random should be set to 1 if W (habitat raster) is NULL. Correcting now...")
    n_random = 1
  }
  
  n_unique = nrow(prev_state) # number of unique paths / velocity-position combinations
  n_r_tot = nrow(prev_state) * n_random
  
  # For calculations that will require shifting and moving of random points
  range_centre_repped_random = apply(range_centre, 2, rep, each = n_random)
  # If this has 1 row, it'll simplify (and simplify = FALSE in "apply" call doesn't fix it) so have to do so manually
  if (!is.matrix(range_centre_repped_random)) range_centre_repped_random = matrix(range_centre_repped_random, 1)
  
  # Get a series of expected values for each velocity and position in each direction
  E_x = as.matrix(prev_state[, c(1, 3), drop = FALSE]) %*% t(H)
  E_y = as.matrix(prev_state[, c(2, 4), drop = FALSE]) %*% t(H)
  
  # Add in random variance based on state covariance
  r_x = apply(E_x, 2, rep, each = n_random) + mvtnorm::rmvnorm(n = n_r_tot, sigma = R)
  r_y = apply(E_y, 2, rep, each = n_random) + mvtnorm::rmvnorm(n = n_r_tot, sigma = R)
  
  res_dt = data.table(t = t_label,
                      path_id = rep(1:n_unique, each = n_random),
                      cpx = r_x[, 1],
                      cpy = r_y[, 1],
                      cvx = r_x[, 2],
                      cvy = r_y[, 2])
  
  if (!is.null(bound)) {
    # Correct any values that land within the boundaries
    locs_sf = st_as_sf(res_dt[, .(cpx, cpy)] * x_scale + range_centre_repped_random, coords = c("cpx", "cpy"), crs = st_crs(bound))
    na_ext_vals = !(1:nrow(locs_sf) %in% st_intersects(bound, locs_sf)[[1]])
    
    # Check if any points fall outside of the boundary. If they're all good, we go ahead without changing anything.
    if (any(na_ext_vals)) {
      n_resim = sum(na_ext_vals)
      # Get nearest points within the boundary for each location that's bad
      nearest_inbounds_pts = st_nearest_points(locs_sf[na_ext_vals, ], bound)
      nearest_pts_coords = st_coordinates(nearest_inbounds_pts)[seq(from = 2, to = n_resim * 2, by = 2), 1:2, drop = FALSE]
      nearest_pts_shifted = (nearest_pts_coords - range_centre_repped_random[na_ext_vals, ]) / x_scale
      # Add in new values
      res_dt$cpx[na_ext_vals] = nearest_pts_shifted[, 1]
      res_dt$cpy[na_ext_vals] = nearest_pts_shifted[, 2]
      
      # Correct velocities using direction adjustment algorithm below
      scaled_velocities = cbind(r_x[na_ext_vals, 2], r_y[na_ext_vals, 2]) * x_scale
      new_velocities = correct_direction_bound_grd(curr_vel = scaled_velocities,
                                                   curr_pos = nearest_pts_coords,
                                                   bound = bound,
                                                   ...) / x_scale
      
      res_dt$cvx[na_ext_vals] = new_velocities[, 1]
      res_dt$cvy[na_ext_vals] = new_velocities[, 2]
    }
  }
  
  if (!is.null(W)) {
    # Do habitat selection
    res_coords = res_dt[, .(cpx, cpy)] * x_scale + range_centre_repped_random
    res_dt$habitat = terra::extract(W, res_coords, ID = FALSE)
    # Correct NA values to 0 if necessary.
    res_dt$habitat[is.na(res_dt$habitat)] = 0
    # If all habitat values are 0 set them to 1, for now. This may have unintended consequences but we'll get to that later.
    res_dt = res_dt[, sum_habitat := sum(habitat), by = path_id][
      , habitat := ifelse(sum_habitat == 0, 1, habitat)]
    # Sample randomly.
    res_dt = res_dt[, .SD[sample(.N, 1, prob = habitat)], by = path_id][, -c("habitat", "sum_habitat")]
  }
  
  # Need to add animal ID and "rep" everything
  res_dt = res_dt[rep(1:.N, each = do_repeat)]
  res_dt$animal_id = rep(1:do_repeat, n_unique)
  setcolorder(res_dt, neworder = c("t", "path_id", "animal_id", "cpx", "cpy", "cvx", "cvy"))
  res_dt
  
}
