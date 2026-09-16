sqpad = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
  # Mandatory RTMB stuff
  "[<-" <- ADoverload("[<-")
  getAll(pars, data_in, warn = FALSE)
  
  # Parameters ----
  
  ## Fixed effects ----
  # log_lambda_q: intensity of perceptibility Poisson process, conditional on availability, for birds at 2-D distance 0 from counter
  # log_alpha_q: distance decay associated with perceptibility process
  # log_D: log of density
  
  # Data ----
  # t_up: time corresponding to upper bin of each detection, measured in time units since beginning of survey
  # r_up: upper distance bin corresponding to each detection
  # count: number of birds with the same time and distance values
  # n_max: maximum pop size for N-mixture model
  # r_fact: 2/3 usually

  # Likelihood function ----
  
  ## transform parameters ----
  r_pred = r_fact * r_up
  lambda_r = exp(log_lambda_q - (exp(log_alpha_q) * r_pred) ^ 2)
  
  p_r = 1 - exp(-lambda_r * t_up)
  
  ## N-mixture model ----
  n_counts = length(r_lo)

  r_up_rep = rep(r_up, each = n_max)
  p_r_rep = rep(p_r, each = n_max)
  count_rep = rep(count, each = n_max)
  n_rep = rep(1:n_max, n_counts)
  
  keep_inds = n_rep >= count_rep
  
  ## Final likelihood step ----
  nll_pop = dpois(n_rep, exp(log_D) * pi * r_up_rep ^ 2)
  nll_det = dbinom(count_rep, n_rep, p_r_rep)
  
  nll_mat = matrix(nll_pop * nll_det * keep_inds, nrow = n_max, ncol = n_counts)
  
  nll_final = -log(apply(nll_mat, 2, sum))
  
  if (debug) return(data.frame(nll_final = nll_final))
  
  sum(nll_final)
  
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
