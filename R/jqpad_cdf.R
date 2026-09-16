# Calculate the cumulative distribution function of a detection using the JQPAD model
#
# r1: distance from counter, lower bound
# r2: distance from counter, upper bound
# t1: time of first detection, lower bound
# t2: time of first detection, upper bound
# rmax: maximum perceptible distance
# tmax: length of survey
# pars: parameters for JQPAD model: [log_lambda_q, log_alpha_q, mu_q, log_sigma_q]
#
# Returns a value of the same length as r and t with probabilities
jqpad_cdf = Vectorize(function(r1,
                               r2,
                               t1,
                               t2,
                               rmax,
                               tmax,
                               pars = numeric(4)) {
  
  lambda_q = exp(pars[1])
  alpha_q = exp(pars[2])
  mu_q = pars[3]
  sigma_q = exp(pars[4])
  
  if (sigma_q == 0) {
    
    f_num = function(rr) rr * (exp(-lambda_q * exp(-alpha_q * rr) * t1) - exp(-lambda_q * exp(-alpha_q * rr) * t2))
    f_den = function(rr) rr * (1 - exp(-lambda_q * exp(-alpha_q * rr) * tmax))
    
    int_num = integrate(f_num, r1, r2)$value
    int_den = integrate(f_den, 0, rmax)$value
    
    return(int_num / int_den)
    
  }
  
  f_num = function(rr) rr * (exp(-lambda_q * exp(-alpha_q * rr) * t1) - exp(-lambda_q * exp(-alpha_q * rr) * t2)) * (pnorm((log(r2) - log(rr) - mu_q) / sigma_q) - pnorm((log(r1) - log(rr) - mu_q) / sigma_q))
  f_den = function(rr) rr * (1 - exp(-lambda_q * exp(-alpha_q * rr) * tmax)) * pnorm((log(rmax) - log(rr) - mu_q) / sigma_q)
  
  int_num = integrate(f_num, 0, rmax)$value
  int_den = integrate(f_den, 0, rmax)$value
  
  int_num / int_den
  
}, vectorize.args = c("r1", "r2", "t1", "t2"))