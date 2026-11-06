# Closed-form ("cf") version of get_offset_la_mc / get_offset_la
#
# Under the joint QPAD model the probability of detecting a bird within duration t and radius d is
#
#   p = integral_0^d 2r/d^2 * (1 - exp(-t * exp(log_lambda - (exp(log_alpha) * r)^2))) dr
#
# Substituting u = r^2 and x = t * lambda * exp(-alpha^2 u) gives an exact solution in terms of the
# entire exponential integral Ein(x) = integral_0^x (1 - exp(-z)) / z dz = E1(x) + log(x) + euler_gamma:
#
#   p = (Ein(x0) - Ein(x1)) / s,   where s = (alpha * d)^2, x0 = t * lambda, x1 = t * lambda * exp(-s)
#
# This is fully vectorized (no Monte Carlo draws, no expand_grid), so it is exact and much faster.
#
# log_lambda: numeric vector, log quasi-availability parameter for each survey
# log_alpha: numeric vector, log quasi-perceptibility parameter for each survey
# duration: survey duration (minutes); scalar or vector
# distance: survey radius (km); scalar or vector
# log: return the log offset?
# nr: ignored; kept so the arguments match get_offset_la_mc
#
# Returns a numeric vector of offsets (log scale by default)
get_offset_la_cf = function(log_lambda,
                            log_alpha,
                            duration = 3,
                            distance = 1,
                            log = TRUE,
                            nr = 100) {

  n_obs = max(length(log_lambda), length(log_alpha), length(duration), length(distance))
  log_lambda = rep_len(log_lambda, n_obs)
  log_alpha = rep_len(log_alpha, n_obs)
  duration = rep_len(duration, n_obs)
  distance = rep_len(distance, n_obs)

  # entire exponential integral; series near zero avoids log(0) / cancellation problems
  ein = function(x) {
    out = x
    small = !is.na(x) & x < 1e-4
    out[small] = x[small] * (1 - x[small] / 4 + x[small] ^ 2 / 18)
    big = !is.na(x) & !small
    out[big] = log(x[big]) - digamma(1)
    # E1(x) < exp(-x) / x, which is negligible for x > 700; calling expint_E1 there is also very slow
    # (it underflows and warns element by element), so only evaluate it where it matters
    mid = big & x <= 700
    out[mid] = out[mid] + suppressWarnings(expint::expint_E1(x[mid]))
    out
  }

  s = (exp(log_alpha) * distance) ^ 2
  x0 = duration * exp(log_lambda)
  x1 = x0 * exp(-s)

  p_detected = (ein(x0) - ein(x1)) / s

  # when s is tiny the difference above loses precision; midpoint approximation is accurate to O(s^2)
  tiny = !is.na(s) & s < 1e-6
  p_detected[tiny] = 1 - exp(-x0[tiny] * exp(-s[tiny] / 2))

  survey_areas = pi * distance ^ 2
  unit_conversion_rates = survey_areas * 100 # convert from km2 to ha - what is the area of the survey region in ha?
  res = p_detected * unit_conversion_rates

  if (log) return(log(res))
  res

}
