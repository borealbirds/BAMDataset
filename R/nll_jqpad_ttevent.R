jqpadtte = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
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
  # r_max: scalar; maximum distance for each bird
  # count: number of birds with the same time and distance values
  # int_weights: Simpson Rule (or other rule) weights for each integral; length of this vector allows us to determine # of terms
  
  # Likelihood function ----
  
  ## transform parameters ----
  lambda_q = exp(log_lambda_q)
  alpha_q = exp(log_alpha_q)

  n_int = length(int_weights) # number of terms to include in each integral
  n_obs = length(t_up)
  
  ## repeat everything ----
  t_up_rep = rep(t_up, n_int)
  t_lo_rep = rep(t_lo, n_int)
  t_max_rep = rep(t_max, n_int)
  r_up_rep = rep(r_up, n_int)
  r_lo_rep = rep(r_lo, n_int)
  r_max_rep = rep(r_max, n_int)
  
  int_weights_rep = rep(int_weights, each = n_obs)
  dr = r_max_rep / (n_int - 1)
  r_int = rep(seq(0, n_int - 1), each = n_obs) / (n_int - 1) * (r_up_rep - r_lo_rep) + r_lo_rep
  
  pnorm_up_fun = pnorm((log(r_up_rep) - log(r_int) - mu_q) / sqrt_2_sigma_q)
  pnorm_lo_fun = pnorm((log(r_lo_rep) - log(r_int) - mu_q) / sqrt_2_sigma_q)
  pnorm_lo_fun[r_lo_rep == 0] = 0
  pnorm_fun = pnorm_up_fun - pnorm_lo_fun
  lambda_r = exp(log_lambda_q - alpha_q * r_int)
  
  ifun_num = int_weights_rep * dr * pnorm_fun * r_int * (exp(-t_lo_rep * lambda_r) - exp(-t_up_rep * lambda_r))
  ifun_den = int_weights_rep * dr * pnorm_fun * r_int * (1 - exp(-t_max_rep * lambda_r))
  
  ifun_num_mat = matrix(ifun_num, n_obs, n_int)
  ifun_den_mat = matrix(ifun_den, n_obs, n_int)
  
  ## summing these effectively approximates the two integrals, which we then take the log of
  nll_num = log(apply(ifun_num_mat, 1, sum))
  nll_den = log(apply(ifun_den_mat, 1, sum))
  
  ## Final likelihood step ----
  nll_pars = -(dnorm(mu_q, log = TRUE) + dnorm(log_sigma_q, mean = -4, log = TRUE)) * count
  nll_final = -nll_num + nll_den + nll_pars
  
  if (debug) return(data.frame(nll_final = nll_final, 
                               nll_num = nll_num,
                               nll_den = nll_den,
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
