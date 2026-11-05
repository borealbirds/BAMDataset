library(tidyverse)

impsamp = function(f,
                   p,
                   q,
                   N = 100,
                   q_domain = c(0, 1)) {
  
  # sample random values from q
  
  rq = runif(N)
  q_cdf = Vectorize(function(pp) integrate(q, 0, pp)$value, vectorize.args = "pp")
  q_quant = Vectorize(function(qq) uniroot(function(t) q_cdf(t) - qq, interval = q_domain)$root, vectorize.args = "qq")
  rr = q_quant(rq)
  
  # calculate f and p and q at each of those values
  
  f_rr = f(rr)
  p_rr = p(rr)
  q_rr = q(rr)
  
  # sum together to approximate integral
  sum(f_rr * p_rr / q_rr) / N
  
}

# int_est = impsamp(f = function(x) x ^ 2.5,
#                   p = function(x) dbeta(x, 0.1, 1),
#                   q = function(x) 2 - 2 * x,
#                   N = 1e4)
# 
# int_exact = integrate(function(x) x ^ 2.5 * dbeta(x, 0.1, 1), lower = 0, upper = 1)$value