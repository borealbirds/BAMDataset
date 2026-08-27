# Generate likelihood profile confidence intervals for a TMB model object. Does not work (yet) for problems with constrained parameters so if you want to use this, need to transform parameters in TMB (e.g., log, logit)
#
# obj: TMB object returned by MakeADFun
# par_ind: integer; which parameter are we profiling? If more than one, gives multiple confidence intervals
# level: number between 0 and 1 representing level of confidence
# parallel: do we run this in parallel?
# h_init: initial step size for algorithm
# h_final: how close do we want to get to the actual lower bound value?
# start_adapt: during each optimization, do we change the initial guess to the optimal non-focal parameters from the previous optimization? This can help speed up each optimization but can also potentially drag the result out of its local minimum and then the algorithm will not work
# max_iter: after how many iterations should we stop looking (even if the algorithm has not converged; prevents infinite loops for messy functions with huge confidence intervals)
# opt_method: optimization method; currently only "nlminb", "L-BFGS-B", "bobyqa" are allowed
# verbose: integer greater than or equal to 0; how much information to provide while the algorithm is running? 0 does nothing; 1 (default) prints beginning of lower and upper for each parameter; 2 prints everything in 1 plus values of focal parameter; 3 prints everything from 2 plus gradient values
#
# Returns a list of data.frames with all the profiles that can easily be converted to a data.frame with confidence intervals using prof_ci_df() below
tmbprofile_ci_manual = function(obj, 
                                par_ind = seq_len(length(obj$par)), 
                                level = 0.95,
                                parallel = FALSE,
                                h_init = 1e-3,
                                h_final = 1e-4,
                                start_adapt = TRUE,
                                max_iter = 50,
                                opt_method = c("nlminb", "BFGS", "bobyqa"),
                                verbose = 1) {
  
  require(foreach)
  require(doParallel)
  
  opt_method = match.arg(opt_method)
  
  fit_fun = function(par, fn, gr, ...) {
    if (opt_method == "nlminb") return(nlminb(par, fn, gr, control = list(eval.max = 1e4, iter.max = 1e4), ...))
    if (opt_method == "BFGS") return(optim(par, fn, gr, method = "BFGS", control = list(maxit = 1e4), ...))
    if (opt_method == "bobyqa") return(nloptr::bobyqa(par, fn, control = list(maxeval = 1e4), ...))
  }
  
  # TO DO: All of this might not work for random effects
  input_obj_mgc = obj$env$tracemgc
  on.exit({obj$env$tracemgc = input_obj_mgc})
  input_obj_par = obj$env$tracepar
  on.exit({obj$env$tracepar = input_obj_par})
  
  obj$env$tracemgc = verbose > 2
  obj$env$tracepar = verbose > 2
  
  if (parallel) {
    `%doint%` = `%dopar%`
  } else {
    `%doint%` = `%do%`
  }
  
  I_PAR = diag(length(obj$par))
  if (is.null(obj$env$random)) {
    OPTIMAL_PARS = obj$env$last.par.best
  } else {
    OPTIMAL_PARS = obj$env$last.par.best[-obj$env$random]
  }
  OPTIMAL_NLL = obj$fn(OPTIMAL_PARS)
  # how much bigger should we let the negative log-likelihood get?
  TARGET_NLL = OPTIMAL_NLL + qchisq(level, 1) / 2
  
  foreach(i = par_ind, .errorhandling = "pass") %doint% {
    
    if (verbose > 0) message("Starting lower profiling for pars[", i, "] = ", names(obj$par[i]))
    
    # Generate new versions of the function and gradient that take 1 less parameter
    new_fn = function(par, pf = OPTIMAL_PARS[i]) {
      obj$fn(par %*% I_PAR[-i, ] + pf %*% I_PAR[i, ])
    }
    
    new_gr = function(par, pf = OPTIMAL_PARS[i]) {
      obj$gr(par %*% I_PAR[-i, ] + pf %*% I_PAR[i, ]) %*% I_PAR[, -i]
    }
    
    # Run binary search algorithm twice to get lower and upper bounds respectively
    
    # lower bound
    result_df_lo = data.frame(iter = 0, h_i = 0, pf_i = OPTIMAL_PARS[i], cur_nll = OPTIMAL_NLL)
    
    h_i = h_init / 2 # we are about to double it; we need something outside the loop though
    pf_i = OPTIMAL_PARS[i]
    opt_other = OPTIMAL_PARS[-i] # use as the initial guess for the optimizer
    cur_nll = OPTIMAL_NLL
    n_iter_lo = 1
    # First part - figure out how large we have to get before we exceed the target. We keep doubling the step size until we get to a value that is larger than the target (or until we run out of iterations)
    while (cur_nll < TARGET_NLL & n_iter_lo < max_iter) {
      h_i = h_i * 2 # double step size
      pf_i = pf_i - h_i # minus because we are doing the lower bound first
      if (!is.finite(new_fn(opt_other, pf = pf_i))) {
        cur_nll = Inf # make the algorithm think we've gone over
        if (verbose > 1) message("Profile value is not finite. Adjusting step size...")
        result_df_lo = rbind(result_df_lo, data.frame(iter = n_iter_lo, h_i = -h_i, pf_i = pf_i, cur_nll = cur_nll))
        n_iter_lo = n_iter_lo + 1
        break # stop increasing the value if we get to an NA somewhere
      }
      if (verbose > 1) message("Iteration ", n_iter_lo, ": par[", i, "] = ", round(pf_i, 5))
      fit = fit_fun(opt_other, new_fn, new_gr, pf = pf_i)
      cur_nll = new_fn(fit$par, pf = pf_i)
      if (!is.finite(cur_nll)) cur_nll = Inf
      if (is.finite(cur_nll) && cur_nll < OPTIMAL_NLL) warning("The profiling algorithm has found a new optimum, better than the input value by ", round(OPTIMAL_NLL - cur_nll, 3), ". Consider running again.")
      if (verbose > 1) message("Profile value: ", round(cur_nll, 2), "; Target value: ", round(TARGET_NLL, 2))
      if (start_adapt & cur_nll < TARGET_NLL) opt_other = fit$par
      result_df_lo = rbind(result_df_lo, data.frame(iter = n_iter_lo, h_i = -h_i, pf_i = pf_i, cur_nll = cur_nll))
      n_iter_lo = n_iter_lo + 1
    }
    # Second part - if we have made it here it implies that we have "overshot" the target and need to start reducing the step size. TO DO: Why is our stopping condition h_final / 2? I know it has to be but need to think about it more...
    while (abs(h_i) > h_final / 2 & n_iter_lo < max_iter) {
      h_i = abs(h_i) / 2 * sign(cur_nll - TARGET_NLL) # reduce the step size by a factor of 2; plus here because if cur_nll > TARGET_NLL we want to try a larger pf_i (we have gone below the lower bound)
      pf_i = pf_i + h_i
      if (verbose > 1) message("Iteration ", n_iter_lo, ": par[", i, "] = ", round(pf_i, 5))
      if (!is.finite(new_fn(opt_other, pf = pf_i))) {
        cur_nll = Inf # Keep as positive because we always want to go back towards the previous value (we know the true confidence bound is on the other side of the NA's)
        if (verbose > 1) message("Profile value is not finite. Adjusting step size...")
        result_df_lo = rbind(result_df_lo, data.frame(iter = n_iter_lo, h_i = h_i, pf_i = pf_i, cur_nll = cur_nll))
        n_iter_lo = n_iter_lo + 1 # not a stopping condition anymore but still relevant for result_df_lo
        next
      }
      fit = fit_fun(opt_other, new_fn, new_gr, pf = pf_i)
      cur_nll = new_fn(fit$par, pf = pf_i)
      if (!is.finite(cur_nll)) cur_nll = Inf
      if (is.finite(cur_nll) && cur_nll < OPTIMAL_NLL) warning("The profiling algorithm has found a new optimum, better than the input value by ", round(OPTIMAL_NLL - cur_nll, 3), ". Consider running again.")
      if (verbose > 1) message("Profile value: ", round(cur_nll, 2), "; Target value: ", round(TARGET_NLL, 2))
      if (start_adapt & cur_nll < TARGET_NLL) opt_other = fit$par
      result_df_lo = rbind(result_df_lo, data.frame(iter = n_iter_lo, h_i = h_i, pf_i = pf_i, cur_nll = cur_nll))
      n_iter_lo = n_iter_lo + 1 # not a stopping condition anymore but still relevant for result_df_lo
    }
    
    # upper bound
    if (verbose > 0) message("Starting upper profiling for pars[", i, "] = ", names(obj$par[i]))
    
    result_df_up = data.frame(iter = 0, h_i = 0, pf_i = OPTIMAL_PARS[i], cur_nll = OPTIMAL_NLL)
    
    h_i = h_init / 2 # we are about to double it; we need something outside the loop though
    pf_i = OPTIMAL_PARS[i]
    opt_other = OPTIMAL_PARS[-i] # use as the initial guess for the optimizer
    cur_nll = OPTIMAL_NLL
    n_iter_up = 1
    # First part - figure out how large we have to get before we exceed the target. We keep doubling the step size until we get to a value that is larger than the target (or until we run out of iterations)
    while (cur_nll < TARGET_NLL & n_iter_up < max_iter) {
      h_i = h_i * 2 # double step size
      pf_i = pf_i + h_i # plus because we are now doing the upper bound
      if (verbose > 1) message("Iteration ", n_iter_up, ": par[", i, "] = ", round(pf_i, 5))
      if (!is.finite(new_fn(opt_other, pf = pf_i))) {
        cur_nll = Inf # make the algorithm think we've gone over
        if (verbose > 1) message("Profile value is not finite. Adjusting step size...")
        result_df_up = rbind(result_df_up, data.frame(iter = n_iter_up, h_i = -h_i, pf_i = pf_i, cur_nll = cur_nll))
        n_iter_up = n_iter_up + 1
        break # stop increasing the value if we get to an NA somewhere
      }
      fit = fit_fun(opt_other, new_fn, new_gr, pf = pf_i)
      cur_nll = new_fn(fit$par, pf = pf_i)
      if (!is.finite(cur_nll)) cur_nll = Inf
      if (is.finite(cur_nll) && cur_nll < OPTIMAL_NLL) warning("The profiling algorithm has found a new optimum, better than the input value by ", round(OPTIMAL_NLL - cur_nll, 3), ". Consider running again.")
      if (verbose > 1) message("Profile value: ", round(cur_nll, 2), "; Target value: ", round(TARGET_NLL, 2))
      if (start_adapt & cur_nll < TARGET_NLL) opt_other = fit$par
      result_df_up = rbind(result_df_up, data.frame(iter = n_iter_up, h_i = h_i, pf_i = pf_i, cur_nll = cur_nll))
      n_iter_up = n_iter_up + 1
    }
    # Second part - if we have made it here it implies that we have "overshot" the target and need to start reducing the step size.
    while (abs(h_i) > h_final / 2 & n_iter_up < max_iter) {
      h_i = abs(h_i) / 2 * sign(cur_nll - TARGET_NLL) # reduce the step size by a factor of 2; minus here because if cur_nll > TARGET_NLL we want to try a smaller pf_i (we have gone above the upper bound)
      pf_i = pf_i - h_i
      if (verbose > 1) message("Iteration ", n_iter_up, ": par[", i, "] = ", round(pf_i, 5))
      if (!is.finite(new_fn(opt_other, pf = pf_i))) {
        cur_nll = Inf # Keep as positive because we always want to go back towards the previous value (we know the true confidence bound is on the other side of the NA's)
        if (verbose > 1) message("Profile value is not finite. Adjusting step size...")
        result_df_up = rbind(result_df_up, data.frame(iter = n_iter_up, h_i = h_i, pf_i = pf_i, cur_nll = cur_nll))
        n_iter_up = n_iter_up + 1 # not a stopping condition anymore but still relevant for result_df_lo
        next
      }
      fit = fit_fun(opt_other, new_fn, new_gr, pf = pf_i)
      cur_nll = new_fn(fit$par, pf = pf_i)
      if (!is.finite(cur_nll)) cur_nll = Inf
      if (is.finite(cur_nll) && cur_nll < OPTIMAL_NLL) warning("The profiling algorithm has found a new optimum, better than the input value by ", round(OPTIMAL_NLL - cur_nll, 3), ". Consider running again.")
      if (verbose > 1) message("Profile value: ", round(cur_nll, 2), "; Target value: ", round(TARGET_NLL, 2))
      if (start_adapt & cur_nll < TARGET_NLL) opt_other = fit$par
      result_df_up = rbind(result_df_up, data.frame(iter = n_iter_up, h_i = -h_i, pf_i = pf_i, cur_nll = cur_nll))
      n_iter_up = n_iter_up + 1 # not a stopping condition anymore but still relevant for result_df_up
    }
    
    rbind(cbind(result_df_lo, side = "lower"), cbind(result_df_up, side = "upper"))
    
  }
  
  # TO DO: Make the list named with the names of parameters? Does this work when we have parameter vectors? Would also have to revise profile_to_table function.
}

# Generate a likelihood slice of a function with respect to a given parameter or set of parameters
#
# obj: TMB object returned by MakeADFun
# par_ind: integer; which parameter are we profiling? If more than one, gives multiple confidence intervals
# bounds: numeric matrix with (# of pars / length(par_ind)) rows and 2 columns; lower and upper bounds for each parameter
# n_slices: integer >= 2; # of different parameter values to try
# parallel: do we run this in parallel?
# par_in: numeric vector; starting parameter values for each slice
# start_adapt: logical; do we change the starting values for each slice to reflect the previous optima?
# verbose: integer greater than or equal to 0; how much information to provide while the algorithm is running? 0 does nothing; 1 (default) prints beginning of lower and upper for each parameter; 2 prints everything in 1 plus values of focal parameter; 3 prints values of other parameters at each optimization, 4 prints everything from 3 plus gradient values
#
# Returns a list of data.frames with all the profiles that can easily be converted to a data.frame with confidence intervals using prof_ci_df() below
tmbprofile_slice_manual = function(obj, 
                                   par_ind = seq_len(length(obj$par)), 
                                   bounds = cbind(rep(0, length(obj$par)), rep(1, length(obj$par))),
                                   n_slices = 10,
                                   par_in = obj$env$last.par.best,
                                   start_adapt = FALSE,
                                   parallel = FALSE,
                                   verbose = 1) {
  
  require(foreach)
  require(doParallel)
  
  # TO DO fix this so you don't have to submit bounds for parameters you don't care about
  if (nrow(bounds) < length(par_ind)) {
    bounds = bounds[rep(1, length(obj$par)), ]
    warning("Not enough bounds have been supplied. Using first row for all parameters")
  }
  
  # TO DO: All of this might not work for random effects
  input_obj_mgc = obj$env$tracemgc
  on.exit({obj$env$tracemgc = input_obj_mgc})
  input_obj_par = obj$env$tracepar
  on.exit({obj$env$tracepar = input_obj_par})
  
  obj$env$tracemgc = verbose > 3
  obj$env$tracepar = verbose > 3
  
  if (parallel) {
    `%doint%` = `%dopar%`
  } else {
    `%doint%` = `%do%`
  }
  
  I_PAR = diag(length(obj$par))
  if (!is.null(obj$env$random)) par_in = par_in[-obj$env$random]
  
  foreach(i = par_ind, .errorhandling = "pass") %doint% {
    
    if (verbose > 0) message("Starting slice for pars[", i, "] = ", names(obj$par[i]))
    
    bound_ind = which(par_ind == i)
    
    # Generate new versions of the function and gradient that take 1 less parameter
    new_fn = function(par, pf = par_in[i]) {
      obj$fn(par %*% I_PAR[-i, ] + pf %*% I_PAR[i, ])
    }
    
    new_gr = function(par, pf = par_in[i]) {
      obj$gr(par %*% I_PAR[-i, ] + pf %*% I_PAR[i, ]) %*% I_PAR[, -i]
    }
    
    par_seq_i = seq(bounds[bound_ind, 1], bounds[bound_ind, 2], length.out = n_slices)
    data.frame(par = par_seq_i,  value = sapply(par_seq_i, function(pf_i) {
      if (verbose > 1) message("---->Optimizing with pars[", i, "] = ", round(pf_i, 3))
      
      if (start_adapt & pf_i != bounds[bound_ind, 1]) {
        opt_init = obj$env$last.par[-i]
      } else {
        opt_init = par_in[-i]
      }
      
      fit = optim(opt_init, new_fn, new_gr, pf = pf_i, method = "BFGS")
      
      if (verbose > 2) message("Optimal non-slice parameters with pars[", i, "] = ", round(pf_i, 3), ": [", paste0(round(fit$par, 3), collapse = ", "), "]")
      
      new_fn(fit$par, pf = pf_i)
    }))
    
  }
}

# Makes a nice confidence interval table from the output of tmbprofile_ci_manual
profile_to_table = function(out) {
  do.call(rbind, lapply(1:length(out), function(ind) {
    # For each parameter
    this_df = out[[ind]]
    
    out = tryCatch({
      this_df_lo = this_df[this_df$side == "lower", ]
      this_df_hi = this_df[this_df$side == "upper", ]
      
      data.frame(par_name = rownames(this_df)[1], est = this_df$pf_i[1], CL = this_df_lo$pf_i[nrow(this_df_lo)], CU = this_df_hi$pf_i[nrow(this_df_hi)])
    }, error = function(e) {
      data.frame(par_name = "NA_FAILED", est = NaN, CL = NaN, CU = NaN)
    })
    
    out
  }))
}
