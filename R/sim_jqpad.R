# n: number of birds simulated
# pars: list of parameters for RQPAD model
# rbins: numeric vector with increasing values and 0 as the first value, distance bins for output data
# tbins: numeric value with increasing values and 0 as the first value, time bins for output data
# return_all: logical; return exact values?
sim_jqpad = function(n,
                     pars = list(),
                     rbins = c(0, 0.1, 1, 10),
                     tbins = c(0, 3, 5, 10),
                     return_all = FALSE) {
  
  tbins_lo = tbins[-length(tbins)]
  tbins_up = tbins[-1]
  rbins_lo = rbins[-length(rbins)]
  rbins_up = rbins[-1]
  
  if (!all(tbins_lo < tbins_up)) stop("Time bins should be ordered from smallest to largest.")
  if (!all(rbins_lo < rbins_up)) stop("Time bins should be ordered from smallest to largest.")
  
  pars_full = fill_empty_list_jqpad(pars)
  getAll(pars_full, warn = FALSE)
  
  r_max = max(rbins)
  t_max = max(tbins)
  
  xydf = data.frame(x = runif(n, min = -r_max, max = r_max), 
                    y = runif(n, min = -r_max, max = r_max)) %>%
    mutate(rr = sqrt(x ^ 2 + y ^ 2),
           logit_dp = rnorm(n(), 0, exp(log_sigma_p)),
           rt = rnorm(n(), sd = exp(log_sigma_q)),
           r_err = exp(log(rr) + mu_q + rt)) %>%
    dplyr::filter(r_err <= r_max) %>%
    mutate(lambda = exp(log_lambda_q - exp(log_alpha_q) * rr + logit_dp),
           tt = rexp(n(), lambda)) %>%
    dplyr::filter(tt <= t_max) %>%
    mutate(t_cut = cut(tt, tbins),
           r_cut = cut(r_err, rbins),
           t_lo = tbins_lo[t_cut],
           t_up = tbins_up[t_cut],
           r_lo = rbins_lo[r_cut],
           r_up = rbins_up[r_cut],
           t_max = as.numeric(max(tbins_up)),
           r_max = as.numeric(max(rbins_up)))
  
  if (!return_all) xydf = xydf %>% dplyr::select(t_lo, t_up, t_max, r_lo, r_up, r_max)
  
  xydf
  
}

# Helper function for above, which adds default parameter values to a list of parameters (that may be missing some)
fill_empty_list_jqpad = function(L,
                                 default = list(log_lambda_q = 0,
                                                log_alpha_q = 0,
                                                mu_q = 0,
                                                log_sigma_q = 0,
                                                log_sigma_p = 0)) {
  
  L_mod = list(log_lambda_q = L$log_lambda_q,
               log_alpha_q = L$log_alpha_q,
               mu_q = L$mu_q,
               log_sigma_q = L$log_sigma_q,
               log_sigma_p = L$log_sigma_p)
  
  for (ind in 1:length(L_mod)) {
    
    if (is.null(L_mod[[ind]])) L_mod[[ind]] = default[[ind]]
    
  }
  
  L_mod
  
}