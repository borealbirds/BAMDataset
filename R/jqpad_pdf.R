# Calculate the probability distribution function of a detection using the JQPAD model
#
# r: distance from counter
# t: time of first detection
# rmax: maximum perceptible distance
# tmax: length of survey
# pars: parameters for JQPAD model: [log_lambda_q, log_alpha_q, mu_q, log_sigma_q]
#
# Returns a value of the same length as r and t with probabilities
jqpad_pdf = Vectorize(function(r,
                               t,
                               rmax,
                               tmax,
                               pars = numeric(4)) {
  
  if (r == 0) return(0)
  
  lambda_q = exp(pars[1])
  alpha_q = exp(pars[2])
  mu_q = pars[3]
  sigma_q = exp(pars[4]) 
  
  scalar = lambda_q / r / sigma_q / sqrt(2 * pi)
  
  f_num = function(rr) rr * exp(-alpha_q * rr - lambda_q * exp(-alpha_q * rr) * t - (log(r) - log(rr) - mu_q)^2 / 2 / sigma_q ^ 2)
  f_den = function(rr) rr * (1 - exp(-tmax * lambda_q * exp(-alpha_q * rr))) * pnorm((log(rmax) - log(rr) - mu_q) / sigma_q)
  
  int_num = integrate(f_num, 0, rmax)$value
  int_den = integrate(f_den, 0, rmax)$value
  
  scalar * int_num / int_den
  
}, vectorize.args = c("r", "t"))