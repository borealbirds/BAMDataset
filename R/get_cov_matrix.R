# Function to get covariance matrix for a LARC fit. Should be sensitive to errors
#
# fit: object returned by fit_jpqpadrc() or equivalent function
#
# Returns a symmetric matrix
get_cov_matrix = function(fit) {
  
  if ("hessian" %notin% names(fit)) stop("the 'fit' object must have an element named 'hessian'; maybe try rerunning fitting function with 'return_hessian = TRUE'?")
  if ("ci_table" %notin% names(fit)) stop("the 'fit' object must have an element named 'ci_table'; maybe try rerunning fitting function with 'return_ci = TRUE'?")
  if (!is.numeric(fit$hessian) || !is.matrix(fit$hessian)) stop("'fit$hessian' must be a numeric matrix")
  if (nrow(fit$hessian) != ncol(fit$hessian)) stop("'fit$hessian' must be a square matrix")
  
  ci95_factor = qnorm(0.975)
  this_cov_mat = try(solve(force_symmetric(this_fit$hessian)))
  
  if (is(this_cov_mat, "try-error") || any(is.na(this_cov_mat))) {
    this_cov_mat = diag(sapply(seq_along(fit$ci_table$est), function(ind) ((fit$ci_table$CU[ind] - fit$ci_table$CL[ind]) / ci95_factor / 2) ^ 2))
  }
  
  this_cov_mat
  
}

force_symmetric = function(M) (M + t(M)) / 2