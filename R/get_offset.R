get_offset = function(lambda_pars,
                      alpha_pars,
                      duration = 3,
                      distance = 1,
                      timeofday = "sunrise",
                      timeofyear = 0,
                      biomass = 0,
                      near_road = 0,
                      aru_01 = 0,
                      log = TRUE) {
  
  timeofday_names = paste0("log_lambda_", timeofday)
  
  log_lambda = lambda_pars[timeofday_names] + lambda_pars["log_lambda_timeofyear"] * timeofyear + lambda_pars["log_lambda_is_aru"] * aru_01
  log_alpha = alpha_pars["log_alpha_"] + alpha_pars["log_alpha_open_closed"] * biomass + alpha_pars["log_alpha_near_road"] * near_road
  
  get_offset_la(log_lambda, log_alpha, duration = duration, distance = distance, log = log)
  
}

get_offset_la = function(log_lambda,
                         log_alpha,
                         duration = 3,
                         distance = 1,
                         log = TRUE) {
  
  # get proportion of birds detected by each survey
  p_detected_inrange = Vectorize(pdetectfun, vectorize.args = c("log_lam", "log_alp", "rmax", "tmax"))(log_lambda, log_alpha, duration, distance, method = "joint")
  
  # estimate unit conversion for each survey
  survey_areas = pi * distance ^ 2
  unit_conversion_rates = survey_areas * 100 # convert from km2 to ha - what is the area of the survey region in ha?
  
  res = p_detected_inrange * unit_conversion_rates
  
  if (log) return(log(res))
  res
  
}

get_offset_la_mc = function(log_lambda,
                            log_alpha,
                            duration = 3,
                            distance = 1,
                            log = TRUE,
                            nr = 100) {
  
  in_dat = data.frame(log_lambda = log_lambda,
                      log_alpha = log_alpha,
                      duration = duration,
                      distance = distance)
  n_obs = nrow(in_dat)
  in_dat$par_id = seq_len(n_obs)
  
  dat_dt = as.data.table(expand_grid(in_dat, r = runif(nr, 0, 1)))[
    , p_detected := 2 * r * (1 - exp(-duration * exp(log_lambda - (exp(log_alpha) * (r * distance)) ^ 2)))][
      , .(p_detected = mean(p_detected)), by = par_id]
  
  survey_areas = pi * in_dat$distance ^ 2
  unit_conversion_rates = survey_areas * 100 # convert from km2 to ha - what is the area of the survey region in ha?
  res = dat_dt$p_detected * unit_conversion_rates
  
  if (log) return(log(res))
  res
  
  
}

# What to do 

## Importance sampling?
get_offset_la_is = function(log_lambda,
                            log_alpha,
                            duration = 3,
                            distance = 1,
                            log = TRUE,
                            nr = 100,
                            is_quant = 0.8) {
  
  in_dat = data.frame(log_lambda = log_lambda,
                      log_alpha = log_alpha,
                      duration = duration,
                      distance = distance)
  n_obs = nrow(in_dat)
  in_dat$par_id = seq_len(n_obs)
  
  dat_dt = as.data.table(expand_grid(in_dat, rq = runif(nr, 0, 1)))[
    , q0 := exp(-(exp(log_alpha) * distance) ^ 2)][
      , r := ifelse(rq < is_quant, rq * 0.8, rq * 0.2 + 0.8)][
        , p_det_num := 2 * r * (1 - exp(-duration * exp(log_lambda - (exp(log_alpha) * r) ^ 2))) / distance ^ 2][
          , p_det_den := ifelse(rq < is_quant, 1.25, 5)][
            , p_detected := p_det_num * p_det_den][
              , .(p_detected = mean(p_detected)), by = par_id]
  
  survey_areas = pi * in_dat$distance ^ 2
  unit_conversion_rates = survey_areas * 100 # convert from km2 to ha - what is the area of the survey region in ha?
  res = dat_dt$p_detected * unit_conversion_rates
  
  if (log) return(log(res))
  res
  
  
}

## Approximate the integral with some function that only requires calculating once (i.e., not at every distance)?

