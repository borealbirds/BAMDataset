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
  
  # get proportion of birds detected by each survey
  p_detected_inrange = Vectorize(pdetectfun, vectorize.args = c("log_lam", "log_alp", "rmax", "tmax"))(log_lambda, log_alpha, duration, distance, method = "joint")
  
  # estimate unit conversion for each survey
  survey_areas = pi * distance ^ 2
  unit_conversion_rates = survey_areas * 100 # convert from km2 to ha - what is the area of the survey region in ha?
  
  res = p_detected_inrange * unit_conversion_rates
  
  if (log) return(log(res))
  res
  
}