# Simulate paths of animal groups, optionally with observation error
#
# n_steps: number of timesteps per path
# n_paths: number of independent paths
# n_individuals: number of individuals in each group
# sigma_v: variability in velocity
# sigma_s: variance associated with individual spacing
# sigma_f: variance associated with foraging spacing
# sigma_o: observation error variance (if equal to 0, no spatial error included)
# beta_v: temporal autocorrelation associated with velocity
# beta_s: temporal aucotorrelation associated with individual spacing
# beta_f: temporal autocorrelation associated with foraging spacing
# omega_v: intensity of range restriction
# W_movement: SpatRaster consisting of movement probabilities; instills habitat selection into the "centroid". Values should all be greater than 0
# W_spacing: SpatRaster consisting of movement probabilities; instills habitat selection into spacing. Values should all be greater than 0
# W_movement_inds: integer vector of length n_steps indicating which index (1-based) of W_movement will be used for each time value.
# W_spacing_inds: integer vector of length n_steps indicating which index (1-based) of W_spacing will be used for each time value.
# n_random: how many random values to simulate at each step, from which we'll pick one based on habitat quality. If W rasters for habitat are NULL then this gets set to 1.
# t_scale: the length of one timestep, in the units associated with the movement parameters
# x_scale: the length of one movement parameter unit, in spatial units
# bound_polygon: geom_sf polygon object representing areas where animals cannot travel. When we extract values, should give a non-NA value for eligible locations.
# x0: numeric matrix with n_paths rows and 2 columns representing initial locations of each group. If NULL, either sets to 0 (in the absence of a spatially explicit region) or samples randomly from the provided region. If it only has one row and n_paths > 1, repeats the same condition for each path.
# v0: numeric matrix with n_paths rows and 2 columns representing initial velocity of each group in each direction. If NULL, sets to 0 for each path. If it only has one row and n_paths > 1, repeats the same condition for each path.
# range_centre: same as x0; where are the animals encouraged to return to? All values are shifted by this value throughout. If NULL, shifts by x0 and that is assumed to be the centre.
# max_range_dist: numeric > 0; how far away are animal(s) allowed to get from home range centre? can be two numbers if it's different in x and y direction
# obs_keep: number between 0 and 1 representing the percentage of observations to be retained
# id_keep: number between 0 and 1 representing the percentage of observations for which individual ID has been retained
# n_cores: integer > 0; if greater than 1, we run in parallel
# do_progress: logical; if TRUE, generates a progress bar indicating how much of the algorithm has completed
# do_plot: logical; if TRUE, plots the result with every timestep. Massively slows down the code but can be useful for debugging or generating animated GIFs of movement paths.
# plot_dir: character representing filepath where plots are to be saved, if they are to be saved. Otherwise, NULL.
# ggplot_extra: ggplot object (e.g., geom_sf(...)) to add on to plot if necessary.
# 
# Returns a data.table with path information
sim_ouf = function(n_steps,
                   n_paths = 1,
                   n_individuals = 1,
                   sigma_v = 1,
                   sigma_s = 0,
                   sigma_f = sigma_s,
                   sigma_o = 0,
                   beta_v = 1,
                   beta_s = 1,
                   beta_f = 1,
                   omega_v = 1,
                   W_movement = NULL,
                   W_spacing = NULL,
                   W_movement_inds = rep(1, n_steps),
                   W_spacing_inds = rep(1, n_steps),
                   n_random = 1,
                   t_scale = 1,
                   x_scale = 1,
                   bound_polygon = NULL,
                   x0 = NULL,
                   v0 = NULL,
                   range_centre = NULL,
                   max_range_dist = Inf,
                   obs_keep = 1,
                   id_keep = 1,
                   n_cores = 1,
                   do_progress = FALSE,
                   do_plot = FALSE,
                   plot_dir = NULL,
                   ggplot_extra = NULL) {
  
  require(tidyverse)
  require(sf)
  require(progress)
  require(data.table)
  require(foreach)
  require(doParallel)
  
  if (n_cores > 1 && do_progress) {
    warning("do_progress cannot be TRUE if n_cores is greater than 1, because progress bars will not look nice when the code is parallelized. Switching off do_progress...")
    do_progress = FALSE
  }
  if (n_cores > 1 && do_plot) {
    warning("do_plot cannot be TRUE if n_cores is greater than 1, because progress bars will not look nice when the code is parallelized. Switching off do_plot")
    do_plot = FALSE
  }
  if (omega_v > beta_v) stop("omega_v must be less than or equal to beta_v. please try again!")
  if (sigma_s > sigma_f) stop("sigma_s must be less than or equal to sigma_f. please try again!")
  
  # Initialize x0 if it has not been initialized already
  if (!is.null(x0) & !is.null(bound_polygon)) {
    # Confirm that initial location is within the polygon and if not, correct it
    x0_sf = st_as_sf(as.data.frame(x0), coords = names(as.data.frame(x0))[1:2], crs = st_crs(bound_polygon))
    x0_out_of_bounds = !(1:n_paths %in% st_intersects(bound_polygon, x0_sf)[[1]])
    
    if (any(x0_out_of_bounds)) {
      n_resim = sum(x0_out_of_bounds)
      # Get nearest points within the boundary for each location that's bad
      nearest_inbounds_pts = st_nearest_points(x0_sf[x0_out_of_bounds, ], bound_polygon)
      nearest_pts_coords = st_coordinates(nearest_inbounds_pts)[seq(from = 2, to = n_resim * 2, by = 2), 1:2, drop = FALSE]
      x0[x0_out_of_bounds, ] = nearest_pts_coords
    }
    
  } else if (is.null(x0) & !is.null(bound_polygon)) {
    # Sample from the appropriate spatial extent, igorning areas within the boundary.
    x0 = st_coordinates(st_sample(bound_polygon, n_paths))
  } else if (is.null(x0)) {
    # Implies that bound_polygon is NULL. In this case, everything is relative so we just set every entry of x0 to 0
    x0 = matrix(0, n_paths, 2)
  }
  if (n_paths > 1 && nrow(x0) == 1) x0 = apply(x0, 2, rep, n_paths)
  
  # Initialize v0 if it has not been initialized already
  if (is.null(v0)) v0 = matrix(0, n_paths, 2)
  if (n_paths > 1 && nrow(v0) == 1) v0 = apply(v0, 2, rep, n_paths)
  
  if (is.null(range_centre)) range_centre = x0 # implies that the animals begin at their range centre wherever that is
  if (nrow(range_centre) == 1 && n_paths > 1) range_centre = apply(range_centre, 2, rep, n_paths)
  x0 = (x0 - range_centre) / x_scale # Shift data
  v0 = v0 / x_scale / t_scale # shift velocity too
  range_centre_repped = apply(range_centre, 2, rep, each = n_individuals) # repeat for each individual now
  if (n_paths == 1) range_centre_repped = matrix(range_centre_repped, 1)
  
  if (length(max_range_dist) > 2) {
    max_range_dist = max_range_dist[1:2]
  } else if (length(max_range_dist) == 1) {
    max_range_dist = rep(max_range_dist, 2)
  }
  check_max_range = any(is.finite(max_range_dist))
  
  # Instantiate some important constants
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
  
  this_H = matrix(nrow = 2, ncol = 2, data = c(H_pp, H_pv, H_vp, H_vv))
  this_R = matrix(nrow = 2, ncol = 2, data = c(V_pp, V_vp, V_vp, V_vv))
  
  E_SS = exp(-t_scale * beta_s / 2)
  E_SF = exp(-t_scale * beta_f / 2)
  SD_S = sigma_s * sqrt(1 - exp(-t_scale * beta_s))
  SD_F = sigma_f * sqrt(1 - exp(-t_scale * beta_f))
  
  # Define this "helper" inside our R function because then we don't need to re-pass all the arguments
  inner_path_gen_func = function(np) {
    out = vector(mode = "list", length = n_steps + 1)
    
    # Make the first timestep based on initial conditions. TO DO: Sample from initial habitat values? Would also have to change format for index vectors
    out[[1]] = data.table(t = 0,
                          path_id = rep(1:np, each = n_individuals),
                          animal_id = rep(1:n_individuals, np),
                          cpx = x0[, 1],
                          cpy = x0[, 2],
                          cvx = v0[, 1],
                          cvy = v0[, 2]) %>% 
      r_spacing(ss_sds = sigma_s,
                sf_sds = sigma_f,
                bound = bound_polygon,
                x_scale = x_scale,
                range_centre = range_centre_repped,
                n_random = n_random,
                W = W_spacing[[W_spacing_inds[1]]])
    
    if(do_progress) pb = progress_bar$new(total = n_steps) 
    
    # Figure out how to name files, if necessary
    if (do_plot && !is.null(plot_dir)) this_digits = paste0("%0", ceiling(log10(n_steps)), "d")
    
    for (tt in 1:n_steps) {
      
      # Select unique values from each path to remove unique individuals within a group, because at this point we only care about the centroid
      new_centroid_locs = r_ouf(prev_state = unique(out[[tt]], by = "path_id")[, .(cpx, cpy, cvx, cvy)],
                                H = this_H,
                                R = this_R,
                                bound = bound_polygon,
                                do_repeat = n_individuals,
                                t_label = tt,
                                x_scale = x_scale,
                                range_centre = range_centre,
                                W = W_movement[[W_movement_inds[tt]]],
                                n_random = n_random)
      
      # Add spacing information
      out[[tt + 1]] = new_centroid_locs %>% 
        r_spacing(ssx_means = out[[tt]]$ssx * E_SS,
                  ssy_means = out[[tt]]$ssy * E_SS,
                  sfx_means = out[[tt]]$sfx * E_SF,
                  sfy_means = out[[tt]]$sfy * E_SF,
                  ss_sds = SD_S,
                  sf_sds = SD_F,
                  bound = bound_polygon,
                  x_scale = x_scale,
                  range_centre = range_centre_repped,
                  W = W_spacing[[W_spacing_inds[[tt]]]],
                  n_random = n_random)
      
      if (check_max_range) {
        new_centroid_range_dx = abs(new_centroid_locs$cpx - range_centre[, 1]) * x_scale
        new_centroid_range_dy = abs(new_centroid_locs$cpy - range_centre[, 2]) * x_scale
        if (all(new_centroid_range_dx > max_range_dist[1] | new_centroid_range_dy > max_range_dist[2])) return(out[seq_len(tt + 1)])
      }
      
      if (do_plot) {
        plot_df = do.call(rbind, out[1:(tt + 1)])[, .(x, y, animal_id)][
          , c("x_scaled", "y_scaled") := .(x * x_scale + range_centre_repped[, 1],
                                           y * x_scale + range_centre_repped[, 2])]
        
        p = ggplot(plot_df)
        
        # For now, fill the in-bounds area with water
        if (!is.null(bound_polygon)) p = p + geom_sf(data = bound_polygon, fill = "#aaddff")
        
        p = p +
          geom_path(aes(x = x_scaled, y = y_scaled, group = animal_id), color = "gray") + 
          geom_point(aes(x = x_scaled, y = y_scaled, color = animal_id)) + 
          ggplot_extra +
          scale_color_viridis_c() +
          ggtitle(paste0("Time: ", as.POSIXct("2025-01-01 12:00:00", tz = "America/Vancouver") + tt * 60 * 60 * t_scale)) + 
          labs(x = "Longitude (degrees)", y = "Latitude (degrees)") +
          custom_ggplot_theme() + 
          theme(legend.position = "none")
        
        if (is.null(plot_dir)) {
          print(p) # Only print plot if it's not writing to a file already (saves a bit of time)
        } else {
          ggsave(paste0(plot_dir, "/t_", sprintf(this_digits, tt), ".png"), plot = p, units = "px", width = 1800, height = 900)
        }
      }
      
      if (do_progress) pb$tick()
      
    }
    
    if (do_progress) pb$terminate()
    
    out
    
  }
  
  if (n_cores > 1) {
    # Run in parallel
    n_paths_per_core = ceiling(n_paths / n_cores)
    if (n_paths_per_core * n_cores != n_paths) {
      warning("The number of paths you requested does not divide evenly into the number of cores you requested, so the algorithm will run ", n_paths_per_core, " paths instead of ", n_paths, ".")
      n_paths = n_paths_per_core * n_cores # correct in case it matters for rest of code
    }
    
    all_steps_list = foreach(ind = 1:n_cores) %dopar% {
      out = inner_path_gen_func(np = n_paths_per_core)
      out %>% mutate(path_id = path_id + n_paths_per_core * (ind - 1))
    }
  } else {
    all_steps_list = inner_path_gen_func(np = n_paths)
  }
  
  # Do observation info after the fact because it is time-independent
  df_final = do.call(rbind, lapply(all_steps_list, r_obs_error, 
                                   sigma_o = sigma_o, 
                                   obs_keep = obs_keep,
                                   id_keep = id_keep, 
                                   bound = bound_polygon,
                                   x_scale = x_scale, 
                                   range_centre = range_centre_repped)) %>%
    # Reorder by time, path ID, whale ID if necessary
    arrange(t, path_id, animal_id)
  
  # Re-shift the data back to the actual coordinates
  df_final = df_final[, c("cvx", "cvy", "cpx", "cpy", "x", "y", "x_obs", "y_obs") := .(cvx * x_scale,
                                                                                       cvy * x_scale,
                                                                                       cpx * x_scale + range_centre_repped[, 1],
                                                                                       cpy * x_scale + range_centre_repped[, 2],
                                                                                       x * x_scale + range_centre_repped[, 1],
                                                                                       y * x_scale + range_centre_repped[, 2],
                                                                                       x_obs * x_scale + range_centre_repped[, 1],
                                                                                       y_obs * x_scale + range_centre_repped[, 2])]
  
}
