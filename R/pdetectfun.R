# Calculates the probability of detecting a bird within a given time and distance interval
#
# log_lam: quasi-availability parameter for model
# log_alp: quasi-perceptibility parameter for model
# tmax: duration of survey
# rmax: radius of spatial area to be surveyed
# method: whether time and distance components of the model should be computed independently (QPAD) or jointly (JQPAD / SQPAD)
# zero_value: what to return when survey duration is set to 0
#
# Returns a numeric value between 0 and 1
pdetectfun = function(log_lam,
                      log_alp,
                      tmax,
                      rmax,
                      method = c("ind", "joint"),
                      zero_value = 0) {
  
  method = method[1]
  
  lam = exp(log_lam)
  alp = exp(log_alp)
  
  if (rmax == 0) return(1 - exp(-lam * tmax))
  if (tmax == 0) return(zero_value)
  
  if (method == "ind") {
    detfun = function(r) 2 * r / rmax ^ 2 * (1 - exp(-lam * tmax)) * exp(-(r * alp) ^ 2)
  } else {
    detfun = function(r) 2 * r / rmax ^ 2 * (1 - exp(-tmax * exp(log_lam - (alp * r) ^ 2)))
  }
  
  tryCatch(integrate(detfun, 0, rmax)$value,
           error = function(e) pracma::quadinf(detfun, 0, rmax)$Q)
  
}
