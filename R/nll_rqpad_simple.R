rqpadexs = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
  # Mandatory RTMB stuff
  "[<-" <- ADoverload("[<-")
  getAll(pars, data_in, warn = FALSE)
  
  # Parameters ----
  
  ## Fixed effects ----
  # log_lambda_q: intensity of perceptibility Poisson process, conditional on availability, for birds at 2-D distance 0 from counter
  # log_alpha_q: distance decay associated with perceptibility process

  ## Random effects ----
  # logit_dp: true (logit-transformed) availability rates, with mean subtracted off
  
  # Data ----
  # tt: time corresponding to each detection, measured in time units since beginning of survey
  # rr: distance corresponding to each detection
  # t_max: last time at which the bird could've been counted during this survey
  
  # Likelihood function ----
  
  ## Transform random effects into their actual forms ----
  lambda = exp(log_lambda_q - exp(log_alpha_q) * rr + logit_dp)

  ## Fixed effects: Probability of observing time-distance relationships ----
  nll_fixed = -dexp(tt, lambda, log = TRUE)
  
  nll_weights = log(-expm1(-lambda * t_max))
  
  ## Random effects: Probability of each random effect taking its value ----
  nll_pp = -dnorm(logit_dp, 0, exp(log_sigma_p), log = TRUE)
  
  ## Final likelihood step ----
  nll_final = nll_fixed + nll_weights + nll_pp
  
  if (debug) return(data.frame(nll_final = nll_final, 
                               nll_fixed = nll_fixed,
                               nll_weights = nll_weights,
                               nll_pp = nll_pp,
                               logit_dp = logit_dp,
                               lambda = lambda,
                               rr = rr,
                               tt = tt))
  
  sum(nll_final)
  
}

# Helper functions for working with and debugging RTMB (may move to another file eventually)
envpar_to_list = function(envpar) {
  
  allnames = unique(names(envpar))
  
  out_list = lapply(allnames, function(nm) {
    
    out = envpar[names(envpar) == nm]
    names(out) = NULL
    out
    
  })
  names(out_list) = allnames
  out_list
  
}
