constrained_gls <- function(y, x, sigma_inv, 
                            len_i, len_t, ...) {
  # Solving constrained GLS (see Eq.8 of Ye et al.)
  # Inputs:
  #       y (dbl vec) vector of response variable
  #       x (sparse mat) design matrix
  #       sigma_inv (dbl mat) inv of cov matrix of y
  #       len_i (int) number of time series (I in the paper)
  #       len_t (int) number of periods in a day (T in the paper)
  # Outputs:
  #       (dbl vec) estimated coefficient
  require(CVXR)
  len_beta <- 7 * len_i * len_t
  x_sigma_x <- Matrix::t(x) %*% sigma_inv %*% x
  y_sigma_x <- t(y) %*% sigma_inv %*% x
  
  c_mat <- map(
    seq(1, len_beta, by = len_t),
    ~c(rep(0, .x - 1), rep(1, 3), rep(0, len_beta - .x - 2))
  )
  
  c_mat <- matrix(
    unlist(c_mat),
    ncol = len_beta,
    byrow = T
  )
  one <- rep(1, nrow(c_mat))
  
  beta <- Variable(len_beta)
  obj <- Minimize(quad_form(beta, x_sigma_x) - 2*y_sigma_x %*% beta)
  constr <- c_mat %*% beta == one
  
  prob <- Problem(
    obj,
    constraints = list(constr)
  )
  out <- solve(prob, ...)
  out <- out$getValue(beta)
  return(out)
}

