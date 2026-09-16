# The probability distribution function of detection times for a given set of parameters
#
# t: time
# log_lam: quasi-availability parameter for model
# log_alp: quasi-perceptibility parameter for model
# tmax: duration of survey
# rmax: radius of spatial area to be surveyed
#
# Returns a numeric
dtdetect = function(t,
                    log_lam,
                    log_alp,
                    tmax,
                    rmax) {
  
  num_fun = function(r) r * exp(log_lam - (exp(log_alp) * r) ^ 2 - t * exp(log_lam - (exp(log_alp) * r) ^ 2))
  den_fun = function(r) r * (1 - exp(-tmax * exp(log_lam - (exp(log_alp) * r) ^ 2)))
  
  integrate(num_fun, 0, rmax)$value / integrate(den_fun, 0, rmax)$value
  
}

ptdetect = function(t,
                    log_lam,
                    log_alp,
                    tmax,
                    rmax) {
  
  library(pracma)
  
  Vectorize(function(tt, log_lam, log_alp, tmax, rmax) {
    num_fun = function(r) r * pexp(tt, exp(log_lam - (exp(log_alp) * r) ^ 2))
    den_fun = function(r) r * pexp(tmax, exp(log_lam - (exp(log_alp) * r) ^ 2))
    
    int_try = integrate(num_fun, 0, rmax)$value / integrate(den_fun, 0, rmax)$value
    if (any(is.na(int_try))) int_try = quadinf(num_fun, 0, rmax)$Q / quadinf(den_fun, 0, rmax)$Q
    
    int_try
  }, vectorize.args = "tt")(t, log_lam, log_alp, tmax, rmax)
  
}