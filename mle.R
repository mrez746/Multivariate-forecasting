mle <- function(alpha, 
                sd_alpha,
                len_i, 
                n1, n2, 
                nu, one, 
                A, omega,
                algo = "NLOPT_GD_STOGO_RAND",
                tol = 1e-06,
                max_eval = 10000) {
  # Function for estimating alpha by MLE
  # Inputs:
  #       alpha_0 (dbl vec)
  #       len_i
  #       n1
  #       n2
  #       nu
  #       one
  #       A
  #       omega
  #       algo (chr) e.g., "NLOPT_GD_STOGO", "NLOPT_GD_STOGO_RAND"
  #       tol (dbl) 
  #       max_eval (int) 
  # Outputs:
  #       (dbl vec) ML estimate of alpha. Order: region ID, wday
  require(nloptr)
  omega_inv <- solve(omega)
  A_omega <- t(A) %*% omega_inv
  A_omega_A <- t(A) %*% omega_inv %*% A
  
  obj_fun <- function(x) {
    # Value of the objective function
    # Inputs:
    #       x (dbl vec)
    # Outputs:
    #       (dbl) value of the objective function at alpha
    a <- split(x, rep(1:7, len_i))
    out <- rep(0, 7)
    
    for (j in 1:7) {
      
      next_j <- if_else(j == 7, 1, j + 1)
      prev_j <- if_else(j == 1, 7, j - 1)
      
      out[j] <- n1[j] * t(a[[j]]) %*% A_omega_A %*% a[[j]] + 
        n2[j] * t(a[[j]]) %*% omega_inv %*% a[[j]] -
        2 * n2[j] * t(a[[prev_j]]) %*% A_omega %*% a[[j]] - 
        2 * t(a[[prev_j]]) %*% A_omega %*% nu[[j]] %*% one[[j]] + 
        2 * t(a[[j]]) %*% omega_inv %*% nu[[j]] %*% one[[j]]
    }
    out <- sum(out)
    return(out)
  }
  
  grad_fun <- function(x) {
    # Calculate the gradient of the obj fun
    # Inputs:
    #       x (dbl vdc)
    # Outputs:
    #       (dbl vec) gradient of obj fun at alpha
    a <- split(x, rep(1:7, len_i))
    out <- vector("list", 7)
    
    for (j in 1:7) {
      
      next_j <- if_else(j == 7, 1, j + 1)
      prev_j <- if_else(j == 1, 7, j - 1)
      
      out[[j]] <- (2 * n1[j] * A_omega_A + 2 * n2[j] * omega_inv) %*% a[[j]] - 
        (2 * n2[j] * t(A_omega)) %*% a[[prev_j]] - 
        (2 * n2[next_j] * A_omega) %*% a[[next_j]] + 
        2 * omega_inv %*% nu[[j]] %*% one[[j]] - 
        2 * A_omega %*% nu[[next_j]] %*% one[[next_j]]
    }
    out <- do.call(cbind, out) %>% t() %>% as.vector()
    return(out)
  }
  
  opts <- list("algorithm" = algo,
               "xtol_rel" = tol,
               "maxeval" = max_eval)
  
  out <- nloptr(
    x0 = alpha, 
    eval_f = obj_fun, 
    eval_grad_f = grad_fun,
    lb = alpha - sd_alpha,
    ub = alpha + sd_alpha,
    opts = opts
  )
  return(out)
}









