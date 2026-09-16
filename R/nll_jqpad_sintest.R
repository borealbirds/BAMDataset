jqpadext = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
  # Mandatory RTMB stuff
  "[<-" <- ADoverload("[<-")
  getAll(pars, data_in, warn = FALSE)
  
  # Parameters ----
  
  ## Fixed effects ----
  # log_lambda_q: intensity of perceptibility Poisson process, conditional on availability, for birds at 2-D distance 0 from counter
  # log_lambda_q_s: time of day variable
  # log_alpha_q: distance decay associated with perceptibility process
  
  # Data ----
  # tt: time corresponding to each detection, measured in time units since beginning of survey
  # rr: distance corresponding to each detection
  # hr: hour of survey
  # t_max: last time at which the bird could've been counted during this survey
  
  # Likelihood function ----
  
  ## Transform random effects into their actual forms ----
  log_lambda = log_lambda_q + log_lambda_q_s * sin(hr * pi / 6) - exp(log_alpha_q) * rr
  lambda = exp(log_lambda)
  
  ## Fixed effects: Probability of observing time-distance relationships ----
  # nll_fixed = lambda * tt - log_lambda #This is more stable but useless if we cannot get the weighs part to be more stable, too
  nll_fixed = -dexp(tt, lambda, log = TRUE)
  
  nll_weights = log(-expm1(-lambda * t_max))
  
  ## Final likelihood step ----
  nll_final = nll_fixed + nll_weights
  
  if (debug) return(data.frame(nll_final = nll_final, 
                               nll_fixed = nll_fixed,
                               nll_weights = nll_weights,
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
