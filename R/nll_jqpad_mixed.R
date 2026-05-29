nll_jqpadmix = function(pars, data_in = TMB_IN_DATA, debug = FALSE) {
  
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
  # int_weights: Simpson Rule (or other rule) weights for each integral; length of this vector allows us to determine # of terms
  # int_nodes: point values (from 0 to 1, needs to be multiplied) for each integral
  # alpha_covs: matrix of covariates for alpha
  # lambda_covs: matrix of covariates for lambda
  
  # Likelihood function ----
  
  ## Transform parameters ----
  n_int = length(int_weights) # number of terms to include in each integral
  
  alpha_q = exp(alpha_covs %*% log_alpha_q)
  log_lambda_q = lambda_covs %*% log_lambda_q
  
  c_time_inds = t_up == t_lo
  c_space_inds = r_up == r_lo
  
  cc_inds = c_time_inds & c_space_inds
  dc_inds = c_time_inds & !c_space_inds
  cd_inds = !c_time_inds & c_space_inds
  dd_inds = !c_time_inds & !c_space_inds
  
  ## Continuous-space continuous-time (CSCT) ----
  
  count_cc = count[cc_inds]
  r_cc = r_lo[cc_inds]
  log_lambda_q_cc = log_lambda_q[cc_inds]
  alpha_q_cc = alpha_q[cc_inds]
  log_lambda_r_cc = log_lambda_q_cc - (alpha_q_cc * r_cc) ^ 2
  
  ll_num_cc = log(r_cc) + log_lambda_r_cc - exp(log_lambda_r_cc) * t_lo[cc_inds]
  
  ## Continuous-space discrete-time (CSDT) ----
  
  count_cd = count[cd_inds]
  r_cd = r_lo[cd_inds]
  log_lambda_q_cd = log_lambda_q[cd_inds]
  alpha_q_cd = alpha_q[cd_inds]
  lambda_r_cd = exp(log_lambda_q_cd - (alpha_q_cd * r_cd) ^ 2)
  
  ll_num_cd = log(r_cd) + log(exp(-lambda_r_cd * t_lo[cd_inds]) - exp(-lambda_r_cd * t_up[cd_inds])) 
  
  ## Discrete-space continuous-time (DSCT) ----
  
  count_dc = count[dc_inds]
  n_obs_dc = sum(dc_inds)
  
  t_rep_dc = rep(t_lo[dc_inds], n_int)
  r_up_rep_dc = rep(r_up[dc_inds], n_int)
  r_lo_rep_dc = rep(r_lo[dc_inds], n_int)
  
  dr_rep_dc = r_up_rep_dc - r_lo_rep_dc
  int_weights_rep_dc = rep(int_weights, each = n_obs_dc) * dr_rep_dc
  r_int_num_dc = rep(int_nodes, each = n_obs_dc) * dr_rep_dc + r_lo_rep_dc
  
  log_lambda_q_rep_dc = rep(log_lambda_q[dc_inds], n_int)
  alpha_q_rep_dc = rep(alpha_q[dc_inds], n_int)
  lambda_r_num_dc = exp(log_lambda_q_rep_dc - (alpha_q_rep_dc * r_int_num_dc) ^ 2)
  dexp_num_fun_dc = lambda_r_num_dc * exp(-lambda_r_num_dc * t_rep_dc)
  
  ifun_num_dc = int_weights_rep_dc * r_int_num_dc * dexp_num_fun_dc
  ifun_num_mat_dc = matrix(ifun_num_dc, n_obs_dc, n_int)
  ll_num_dc = log(apply(ifun_num_mat_dc, 1, sum))
  
  ## Discrete-space discrete-time (DSDT) ----
  
  count_dd = count[dd_inds]
  n_obs_dd = sum(dd_inds)
  
  ## repeat everything ----
  t_up_rep_dd = rep(t_up[dd_inds], n_int)
  t_lo_rep_dd = rep(t_lo[dd_inds], n_int)
  r_up_rep_dd = rep(r_up[dd_inds], n_int)
  r_lo_rep_dd = rep(r_lo[dd_inds], n_int)
  
  dr_rep_dd = r_up_rep_dd - r_lo_rep_dd
  int_weights_rep_dd = rep(int_weights, each = n_obs_dd) * dr_rep_dd
  r_int_num_dd = rep(int_nodes, each = n_obs_dd) * dr_rep_dd + r_lo_rep_dd
  
  log_lambda_q_rep_dd = rep(log_lambda_q[dd_inds], n_int)
  alpha_q_rep_dd = rep(alpha_q[dd_inds], n_int)
  lambda_r_num_dd = exp(log_lambda_q_rep_dd - (alpha_q_rep_dd * r_int_num_dd) ^ 2)
  pexp_num_fun_dd = exp(-lambda_r_num_dd * t_lo_rep_dd) - exp(-lambda_r_num_dd * t_up_rep_dd)
  
  ifun_num_dd = int_weights_rep_dd * r_int_num_dd * pexp_num_fun_dd
  ifun_num_mat_dd = matrix(ifun_num_dd, n_obs_dd, n_int)
  ll_num_dd = log(apply(ifun_num_mat_dd, 1, sum))
  
  ## Denominator of likelihood function (same for all models) ----
  
  n_obs_den = length(cc_inds)
  
  t_max_rep = rep(t_max, n_int)
  r_mxs_rep = rep(r_mxs, n_int)
  
  int_weights_rep_den = rep(int_weights, each = n_obs_den) * r_mxs_rep
  r_int_den = rep(int_nodes, each = n_obs_den) * r_mxs_rep
  
  log_lambda_q_rep_den = rep(log_lambda_q, n_int)
  alpha_q_rep_den = rep(alpha_q, n_int)
  lambda_r_den = exp(log_lambda_q_rep_den - (alpha_q_rep_den * r_int_den) ^ 2)
  pexp_den_fun = -expm1(-lambda_r_den * t_max_rep)
  # pexp_den_fun = 1 - exp(-lambda_r_den * t_max_rep)
  
  ifun_den = int_weights_rep_den * r_int_den * pexp_den_fun
  ifun_den_mat = matrix(ifun_den, n_obs_den, n_int)
  ll_denom = log(apply(ifun_den_mat, 1, sum))
  
  ## Bring it all together! ----
  nll_final = -(sum(ll_num_cc * count_cc) + sum(ll_num_cd * count_cd) + sum(ll_num_dc * count_dc) + sum(ll_num_dd * count_dd) - sum(ll_denom * count))
  
  if (debug) {
    return(list(nll_final = nll_final,
                cc_inds = cc_inds,
                cd_inds = cd_inds,
                dc_inds = dc_inds,
                dd_inds = dd_inds,
                ll_num_cc = ll_num_cc,
                ll_num_cd = ll_num_cd,
                ll_num_dc = ll_num_dc,
                ll_num_dd = ll_num_dd,
                ll_denom = ll_denom,
                count = count))
  }
  
  nll_final
  
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

# Makes RTMB model object for this likelihood function. Note that formulas do not currently accept random effects or interaction terms, only "1", "0" (no intercept), or column names of "data_in". If "1" or "0" are not specified assumes intercept is present
#
# data_in: data.frame, assumed to already have the names t_lo, t_up, t_max, r_lo, r_up, r_mxs, count
# formula_lambda: formula object describing which covariates should go into lambda
# formula_alpha: formula object describing which covariates should go into alpha
# n_integral: integer > 0; for integral approximation
# fct_integral: numeric between 0 and 1 (exclusive); for integral approximation
# fit: logical; if TRUE, fits the model as well as returning the RTMB model object
# return_data: logical; if TRUE, returns the data object used for RTMB
# ...: additional arguments to mle() (only relevant if fit is TRUE)
#
# Returns an RTMB model object if fit = FALSE as well as a model fit if fit = TRUE
fit_jqpadmix = function(data_in,
                        formula_lambda = ~ 1,
                        formula_alpha = ~ 1,
                        n_integral = 150,
                        fct_integral = 0.96,
                        fit = TRUE,
                        return_data = FALSE,
                        ...) {
  
  terms_lambda = terms(formula_lambda)
  terms_alpha = terms(formula_alpha)
  
  # extract covariate names
  lambda_cov_names = attr(terms_lambda, "term.labels")
  alpha_cov_names = attr(terms_alpha, "term.labels")
  
  # get covariate values from data
  this_lambda_covs = as.matrix(data_in[, lambda_cov_names])
  this_alpha_covs = as.matrix(data_in[, alpha_cov_names])
  
  # remove covariates that have the same value everywhere
  lambda_cov_sds = apply(this_lambda_covs, 2, sd)
  alpha_cov_sds = apply(this_alpha_covs, 2, sd)
  
  this_lambda_covs = this_lambda_covs[, lambda_cov_sds > 0, drop = FALSE]
  this_alpha_covs = this_alpha_covs[, alpha_cov_sds > 0, drop = FALSE]
  
  # add intercept if necessary
  lambda_intercept = attr(terms_lambda, "intercept")
  alpha_intercept = attr(terms_alpha, "intercept")
  
  if (lambda_intercept > 0) this_lambda_covs = cbind(1, this_lambda_covs)
  if (alpha_intercept > 0) this_alpha_covs = cbind(1, this_alpha_covs)
  
  this_init_par = list(log_lambda_q = numeric(ncol(this_lambda_covs)),
                       log_alpha_q = numeric(ncol(this_alpha_covs)))
  
  integral_info = qfun(n_integral, fct_integral)
  
  rtmb_data_in = data_in %>%
    dplyr::select(r_lo, r_up, r_mxs, t_lo, t_up, t_max, count) %>%
    as.list %>%
    c(list(int_weights = integral_info$weights,
           int_nodes = integral_info$nodes,
           lambda_covs = this_lambda_covs,
           alpha_covs = this_alpha_covs))
  this_nll = function(pars) nll_jqpadmix(pars, data_in = rtmb_data_in)
  
  obj = RTMB::MakeADFun(this_nll, parameters = this_init_par)
  
  if (!fit && !return_data) return(obj)
  if (!fit) return(list(obj = obj, dat = rtmb_data_in))
  if (!return_data) return(list(obj = obj, fit = mle(obj, ...)))

  list(obj = obj,
       fit = mle(obj, ...),
       dat = rtmb_data_in)
  
}

qfun = function(n_pts = 10, fct = 0.5) {
  
  bound_up = fct ^ seq(0, n_pts - 1)
  bound_lo = c(bound_up[-1], 0)
  
  bound_mid = (bound_lo + bound_up) / 2
  bound_diff = bound_up - bound_lo
  
  list(nodes = bound_mid,
       weights = bound_diff)
  
}

int_q = function(l = 0, a = 0, TM = 3, RM = 1, ...) {
  
  ff = function(x) -expm1(-TM * exp(l - (exp(a) * x) ^ 2)) * x
  
  q_info = qfun(...)
  
  sum(ff(q_info$nodes * RM) * q_info$weights * RM)
  
}

int_q_ex = function(l = 0, a = 0, TM = 3, RM = 1) {
  
  ff = function(x) -expm1(-TM * exp(l - (exp(a) * x) ^ 2)) * x
  
  integrate(ff, lower = 0, upper = RM)$value
  
}


err_q = function(l = 0, a = 0, TM = 3, RM = 1, rel = FALSE, ...) {
  
  val_q = int_q(l, a, TM, RM, ...)
  val_ex = int_q_ex(l, a, TM, RM)
  abs_err = val_q - val_ex
  
  if (!rel) return(abs_err)
  
  abs_err / val_ex
  
}

# vals = expand.grid(n = 150, alp = seq(2, 6, 0.05), fct = seq(0.9, 0.999, 0.001)) %>%
#   mutate(err_rel = mapply(err_q, n_pts = n, fct = fct, a = alp, l = 6, rel = TRUE),
#          abs_err_rel = abs(err_rel) * 100,
#          log_abs_err_rel = log10(abs_err_rel))
# 
# vals_rast = vals %>%
#   dplyr::select(x = alp, y = fct, z = log_abs_err_rel) %>%
#   rast()
# 
# plot(vals_rast, asp = diff(range(vals$alp)) / diff(range(vals$fct)))
