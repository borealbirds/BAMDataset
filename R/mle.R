# Optimization function that attempts to reduce the likelihood of ending up in local optima by re-running the profiling algorithm
#
# obj: RTMB model object
# method: optimization method; currently only "nlminb", "L-BFGS-B", "bobyqa" are allowed
# verbose: extent of print output provided. If 0, nothing at all, if 1, minimal output, if 2, prints parameters for every function call, also see tmbprofile_ci_manual which receives this argument
# profile_max: integer > 0; number of times to profile before stopping
# profile_improve_stop: numeric >= 0; how much better does the new likelihood value have to be to keep going?
# eval_max: integer > 0; number of function evaluations allowed during each optimization, passed to internal model fitting functions
# iter_max: integer > 0; number of algorithm iterations allowed during each optimization (NOT for the overall fitting algorithm)
# return_ci: logical; if TRUE, returns a table with confidence intervals
#
# Returns a model fit from nlminb() or some other function (see "method")
mle = function(obj,
               init_par = obj$par,
               method = c("nlminb", "BFGS", "bobyqa"),
               verbose = 0,
               profile_max = 10,
               profile_improve_stop = 0,
               eval_max = 1e5,
               iter_max = 1e5,
               return_ci = FALSE,
               ...) {
  
  # Within this function we will change the amount of print output produced by the TMB model object according to the user supplied verbose argument. To save the original values (thus not altering the user-facing behavior of the object outside of this function), we use an on.exit() call to set them back once the function is done running.
  input_obj_mgc = obj$env$tracemgc
  on.exit({obj$env$tracemgc = input_obj_mgc})
  input_obj_par = obj$env$tracepar
  on.exit({obj$env$tracepar = input_obj_par})
  
  method = match.arg(method)
  
  profile_improve_stop = pmax(profile_improve_stop, 0)
  
  fit_fun = function(par) {
    if (method == "nlminb") return(nlminb(par, obj$fn, obj$gr, control = list(eval.max = eval_max, iter.max = iter_max)))
    if (method == "BFGS") return(optim(par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = iter_max)))
    if (method == "bobyqa") return(nloptr::bobyqa(par, obj$fn, control = list(maxeval = eval_max)))
  }
  
  # Adjust print output preferences according to "verbose"
  obj$env$tracemgc = (verbose > 0)
  obj$env$tracepar = (verbose > 1)
  
  # First attempt at getting best parameters. Doing this will save a (hopefully informed) value under obj$env$last.par.best which tmbprofile_ci_manual() uses by default.
  first_fit = fit_fun(obj$par)
  
  # Calculate profile likelihood confidence intervals. This method surveys the parameter space in a way that is a bit more broad than an optimizer. Often, in the process of profiling the algorithm uncovers a new, better local minimum (potentially the global minimum).
  run_profile = TRUE
  n_profile = 0
  while(run_profile & n_profile < profile_max) {
    
    n_profile = n_profile + 1
    best_value = obj$fn(obj$env$last.par.best)
    if (verbose > 0) message("beginning profile number ", n_profile)
    pr = tmbprofile_ci_manual(obj, verbose = verbose, opt_method = method)
    
    # Compare "best_value" (the old MLE) to whatever the new best value is according to the function. If it's better, we run the profile again to see if we can further determine new optima.
    new_best_value = obj$fn(obj$env$last.par.best)
    run_profile = new_best_value < (best_value - profile_improve_stop)
    
  }
  
  # Fit the model, now starting at the best parameter value, just in case some neighboring parameter combination is a little better. Return this as the function output.
  out = fit_fun(obj$env$last.par.best)
  
  # If we are to return a confidence interval table, generate that table here
  if (return_ci) {
    
    # if possible, we don't want to rerun the profile because it can be slow. if the function value at the optimum is similar to the best value identified by the final optimization, we just use the last profile from the algorithm above. if there is no profile, we run one.
    rerun_profile = is(try(pr), "try-error") || (obj$fn(obj$env$last.par.best) < (pr[[1]]$cur_nll[1] - profile_improve_stop))
    if (rerun_profile) {
      pr = tmbprofile_ci_manual(obj, verbose = verbose, opt_method = method)
    }
    
    # convert to table and replace the "est" column with the actual MLE in case they differ
    pr_table = profile_to_table(pr)
    pr_table$est = obj$env$last.par.best
    
    out = list(fit = out, ci_table = pr_table)
    
  }
  
  out
  
}