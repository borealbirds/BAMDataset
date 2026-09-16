rqpad = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
  # Mandatory RTMB stuff
  "[<-" <- ADoverload("[<-")
  getAll(pars, data_in, warn = FALSE)
  
  # Parameters ----
  
  ## Fixed effects ----
  # log_lambda_q: intensity of perceptibility Poisson process, conditional on availability, for birds at 2-D distance 0 from counter
  # log_alpha_q: distance decay associated with perceptibility process
  # mu_q: bias in log-distance estimated by counter
  # log_sigma_q: variance (around potentially biased mean) in log-distance estimated by counter
  # logit_mu_p: mean availability proportion at an infinite time horizon
  # log_gamma_p: time decay associated with variance of availability proportions within samples
  # log_sigma_p: baseline variance of availability proportions within samples
  
  ## Random effects ----
  # log_dr: difference between "true" distances identified by counter and actual "true" distances, in log-space
  # logit_rt_transf: "true" distances "identified" by the counter
  # logit_tt_transf: true (logit-transformed within user-provided time intervals) times of first detection by the counter
  # logit_dp: true (logit-transformed) availability rates, with mean subtracted off
  
  # Data ----
  # t_lo: lower bound of time bin corresponding to each detection, measured in time units since beginning of survey
  # t_up: upper bound of time bin corresponding to each detection, measured in time units since beginning of survey
  # r_lo: lower bound of distance bin corresponding to each detection
  # r_up: upper bound of distance bin corresponding to each detection
  
  # Likelihood function ----
  
  ## Transform random effects into their actual forms ----
  rt = pnorm(logit_rt_transf) * (r_up - r_lo) + r_lo
  rr = rt * exp(mu_q + log_dr)
  tt = pnorm(logit_tt_transf) * (t_up - t_lo) + t_lo
  pp = 1 / (1 + exp(-(logit_dp + logit_mu_p)))

  ## Fixed effects: Probability of observing time-distance relationships ----
  nll_fixed = -dexp(tt, exp(log_lambda_q - exp(log_alpha_q) * rr) / pp, log = TRUE)
  
  ## Random effects: Probability of each random effect taking its value ----
  nll_rr = -dnorm(log_dr, 0, exp(log_sigma_q), log = TRUE)
  nll_rt = -dnorm(logit_rt_transf, log = TRUE)
  nll_tt = -dnorm(logit_tt_transf, log = TRUE)
  nll_pp = -dnorm(logit_dp, 0, sqrt(exp(2 * log_sigma_p) * tt ^ (-exp(log_gamma_p))), log = TRUE)
  
  ## Final likelihood step ----
  nll_final = nll_fixed + nll_rr + nll_rt + nll_tt + nll_pp
  
  if (debug) return(data.frame(nll_final = nll_final, 
                               nll_fixed = nll_fixed,
                               nll_rr = nll_rr,
                               nll_rt = nll_rt,
                               nll_tt = nll_tt,
                               nll_pp = nll_pp,
                               log_dr = log_dr,
                               logit_rt_trasnf = logit_rt_transf,
                               logit_tt_transf = logit_tt_transf,
                               logit_dp = logit_dp,
                               rt = rt,
                               rr = rr,
                               tt = tt,
                               pp = pp))
  
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
