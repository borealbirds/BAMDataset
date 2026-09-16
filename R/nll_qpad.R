qpad = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
  # Mandatory RTMB stuff
  "[<-" <- ADoverload("[<-")
  getAll(pars, data_in, warn = FALSE)
  
  # Parameters ----
  
  ## Fixed effects ----
  # log_lambda_q: intensity of perceptibility Poisson process, conditional on availability, for birds at 2-D distance 0 from counter
  # log_alpha_q: distance decay associated with perceptibility process
  
  # Data ----
  # t_up: time corresponding to upper bin of each detection, measured in time units since beginning of survey
  # t_lo: time corresponding to lower bin of each detection, measured in time units since beginning of survey
  # r_up: upper distance bin corresponding to each detection
  # r_lo: lower distance bin corresponding to each detection
  # t_max: last time at which the bird could've been counted during this survey
  # r_mxs: maximum radius tracked by each survey, NOT r_max which is the maximum perceptual range of the counter
  # count: number of birds with the same time and distance values

  # Likelihood function ----
  
  ## transform parameters ----
  alpha_q = exp(log_alpha_q)
  lambda_q = exp(log_lambda_q)
  
  ## availability ----
  nll_q = (exp(-t_lo * lambda_q) - exp(-t_up * lambda_q)) / (1 - exp(-t_max * lambda_q))
  
  ## detectability ----
  nll_p = (exp(-(r_lo * alpha_q)^2) - exp(-(r_up * alpha_q)^2)) / (1 - exp(-(r_mxs * alpha_q)^2))

  ## Final likelihood step ----
  nll_final = -log(nll_q) - log(nll_p)
  
  if (debug) return(data.frame(nll_final = nll_final,
                               nll_q = nll_q,
                               nll_p = nll_p,
                               count = count))
  
  sum(nll_final * count)
  
}

# Helper function to generate necessary terms for Simpson's rule. Does not account (or ask) for interval width so that should be factored in later
get_weights_simpson = function(nt) {
  
  if ((nt %% 2) != 1 || (nt < 5)) stop("nt should be an odd integer greater than or equal to 5. Please try again!")
  
  rep_ind = (nt - 3) / 2
  sqnc = c(1, rep(c(4, 2), rep_ind), 4, 1)
  
  sqnc / 3
  
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
