multivariate_forecasting <- function(df, 
                                     horizon,
                                     max_iter = 200, 
                                     algo = "NLOPT_GD_STOGO_RAND",
                                     verbose = TRUE) {
  # Description:
  #       Code for implementation of the multivariate model in
  #       Ye, Han, James Luedtke, and Haipeng Shen. 2019. 
  #       “Call Center Arrivals: When to Jointly Forecast Multiple Streams?” 
  #       Production and Operations Management 28 (1): 27–42.
  # Inputs:
  #       df (tibble) w/ following columns
  #         $stream (int) stream ID index 
  #         $call_volume (dbl) number of call volume
  #         $wd (int) weekday index
  #         $d (int) day index, starting from 1
  #         $t (int) intraday time interval index, starting from 1
  #       horizon (int) forecast horizon, i.e., # forecasts to make 
  #       max_iter (int) maximum iteration for iterative estimation procedures
  #       algo (chr) name of algorithm to use in nloptr, for example
  #         "NLOPT_LD_LBFGS"
  #         "NLOPT_GD_STOGO"
  #         "NLOPT_GD_STOGO_RAND"
  #       verbose (lgl) print out iteration process?
  # Outputs:
  #       (list) w/ following components
  #         $df_pred (tibble) w/ following columns
  #             $stream (int) stream ID index
  #             $h (int) forecast horizon
  #             $pred (int) call volume forecast
  #         $step1_converge (lgl) did procedure of Sec 3.2.1 converge?
  #         $step2_converge (lgl) did procedure of Sec 3.2.2 converge?
  #         $params (list) estimated parameters:
  #             $alpha (dbl vec) order: region ID, wday
  #             $f_vec (dbl vec) order: region ID, wday, time interval
  #             $u_vec (dbl vec) order: region ID, day ID
  #             $A (mat) 
  #             $omega (mat) 
  #             $sigma (mat) 
  # Required packages:
  #       MTS package
  #       lubridate package
  #       tidyverse package
  #       Matrix package
  #       CVXR package
  #       nloptr package
  # Dependencies:
  #       initialization_functions.R
  #       mle.R
  #       constrained_gls.R
  # Note:
  #       The code was written on 5 Sep 2023. In the case that you encounter
  #       issues due to package updates, we recommend using the groundhog
  #       package to access the version of packages on 5 Sep 2023.
  require(MTS)
  require(lubridate)
  require(tidyverse)
  
  
  df <- df %>% 
    arrange(stream, d, t) %>% 
    group_by(stream) %>% 
    nest() %>% 
    mutate(
      call_volume = map(data, "call_volume"), 
      wd = map(data, "wd"), 
      d = map(data, "d"), 
      t = map(data, "t")
    ) %>% 
    select(-data) %>% 
    ungroup()
  stopifnot(df$d[[1]][1] == 1)
  stopifnot(df$t[[1]][1] == 1)
  hours <- df$t[[1]] %>% unique() %>% sort()
  stream <- df$stream
  len_i <- nrow(df)
  len_t <- length(hours)
  len_d <- df$d[[1]] %>% max()
  
  df <- df %>% 
    mutate(
      x = map(call_volume, ~sqrt(.x + 1/4))
    )
  
  # Estimating u, f, and sigma (Sec. 3.2.1) ####
  # Step 1 (initializing)
  # 1-a) Initialize the f_{wd, t}^{(i)} parameter
  list_fhat <- pmap(
    df %>% 
      select(x, wd, t) %>% 
      mutate(hours = list(hours)),
    initialize_fhat
  )
  
  df <- df %>% 
    mutate(
      fhat = list_fhat
    )
  
  # 1-b) Initialize the u_d^{(i)} parameters
  list_uhat <- pmap(
    df %>% 
      select(x, fhat, wd, d, t) %>% 
      mutate(hours = list(hours)),
    initialize_ud
  )
  
  df <- df %>% 
    mutate(
      u_d = map(list_uhat, "u_d"),
      res = map(list_uhat, "res")
    )
  
  # 1-c) Initialize the covariance matrix, Sigma
  mat_res <- matrix(
    df$res %>% unlist(), 
    ncol = nrow(df), 
    byrow = F
  )
  sigma <- var(mat_res)
  
  # Step 2 (iterative updating) 
  u_vec <- df$u_d %>% unlist()
  x_vec <- df$x %>% unlist()
  len_x <- length(x_vec)
  d <- df$d %>% unlist()
  wd <- df$wd %>% unlist()
  t <- df$t %>% unlist()
  i <- rep(1:nrow(df), each = length(df$x[[1]]))
  
  # Initialize the f vector
  f_vec <- rep(0, len_i * 7 * len_t)
  
  # Construct a matrix, u_mat, for solving Eq (8) of Ye et al.
  indx_c1 <- 7 * (i - 1) * len_t + 
    (wd - 1) * len_t + t
  
  u_mat <- Matrix::sparseMatrix(
    i = 1:len_x, 
    j = indx_c1, 
    x = u_vec
  )
  
  for (run in 1:max_iter) {
    
    if (verbose) {
      print(paste0("Estimating f, u, and sigma -- Iteration ", run))
    }
    # Construct the cov matrix for solving the GLS
    indx_r <- map(
      1:len_x,
      ~rep(.x, len_i)
    )
    indx_r <- unlist(indx_r)
    
    indx_c2 <- map2(
      d,
      t,
      ~(.x - 1) * len_t + .y + 0:(len_i - 1) * len_d * len_t
    )
    indx_c2 <- unlist(indx_c2)
    
    values <- map(
      i,
      ~sigma[.x, ]
    )
    values <- unlist(values)
    
    cov_mat <- Matrix::sparseMatrix(
      i = indx_r,
      j = indx_c2,
      x = values
    )
    cov_mat_inv <- Matrix::solve(cov_mat)
    
    # 2-a) Estimate the f vector
    f_old <- f_vec
    f_vec <- constrained_gls(
      y = x_vec, 
      x = u_mat, 
      sigma_inv = cov_mat_inv, 
      len_i = len_i, 
      len_t = len_t
    )
    f_vec <- as.vector(f_vec)
    
    # 2-b) Update the sigma matrix
    res <- x_vec - u_mat %*% f_vec
    
    mat_res <- matrix(
      res, 
      ncol = len_i, 
      byrow = F
    )
    sigma <- var(mat_res)
    
    # Construct a matrix, f_mat, to solve Step (c) pg 7 of Ye et al.
    indx_c3 <- (i - 1) * len_d + d
    
    indx_f <- (i - 1) * 7 * len_t + 
      (wd - 1) * len_t + t
    
    f_mat <- Matrix::sparseMatrix(
      i = 1:len_x, 
      j = indx_c3, 
      x = f_vec[indx_f]
    )
    
    values <- map(
      i,
      ~sigma[.x, ]
    )
    values <- unlist(values)
    
    cov_mat <- Matrix::sparseMatrix(
      i = indx_r,
      j = indx_c2,
      x = values
    )
    
    # 2-c) Estimate the U vector
    u_vec <- solve(Matrix::t(f_mat) %*% cov_mat %*% f_mat) %*% 
      (Matrix::t(f_mat) %*% cov_mat %*% x_vec)
    u_vec <- as.vector(u_vec)
    
    # Update the sigma matrix
    res <- x_vec - f_mat %*% u_vec
    
    mat_res <- matrix(
      res, 
      ncol = len_i, 
      byrow = F
    )
    sigma <- var(mat_res)
    
    # Check stopping criteria
    df_criteria <- tibble(
      f_old = f_old,
      f_vec = f_vec,
      wd = rep(1:7, len_i * len_t),
      i = rep(1:len_i, each = 7 * len_t)
    )
    
    df_criteria <- df_criteria %>% 
      mutate(
        e = f_old - f_vec
      ) %>% 
      group_by(i, wd) %>% 
      summarise(
        e = sqrt(sum(e^2) / n()),
        .groups = "drop"
      )
    criteria_reached <- all(df_criteria$e < 10e-8)
    if(criteria_reached) break
    
    # Redefine the u_mat matrix with updated estimates of the U vector
    indx_u <- (i - 1) * len_d + d
    
    u_mat <- Matrix::sparseMatrix(
      i = 1:len_x, 
      j = indx_c1, 
      x = u_vec[indx_u]
    )
  }
  step1_converge <- run < max_iter
  
  # Estimating alpha, A, and omega (Sec. 3.2.2) ####
  indx <- (i - 1) * len_d + d
  ud <- u_vec[indx]
  ud <- ud[t == 1]
  wd <- wd[t == 1]
  i <- i[t == 1]
  
  # Initialize alpha
  params_0 <- initialize_alpha(
    len_i = len_i, 
    i = i, 
    wd = wd, 
    ud = ud
  )
  wd_2 <- wd[i == 1]
  u_mat <- params_0$u_mat
  alpha_mat <- params_0$alpha_mat
  alpha <- params_0$alpha
  A <- matrix(0, nrow = len_i, ncol = len_i)
  omega <- matrix(0, nrow = len_i, ncol = len_i)
  
  for (run in 1:max_iter) {
    
    if (verbose) {
      print(paste0("Estimating alpha, A, and omega -- Iteration ", run))
    }
    # Estimate A and omega
    mod <- VAR(
      x = t(u_mat - alpha_mat), 
      p = 1, 
      include.mean = F,
      output = FALSE
    )
    A_old <- A
    A <- mod$Phi
    omega_old <- omega
    omega <- mod$Sigma
    mu_mat <- A %*% u_mat[, -nrow(u_mat)] - u_mat[, -1]
    nu_list <- map(1:7, ~mu_mat[, wd_2[-1] == .x])
    
    # Estimate alpha
    
    mle_out <- mle(
      alpha = alpha, 
      sd_alpha = params_0$sd_alpha,
      n1 = params_0$n1_vec, 
      n2 = params_0$n2_vec, 
      nu = nu_list, 
      one = params_0$one_list, 
      len_i = len_i, 
      A = A, 
      omega = omega, 
      tol = 1.0e-6,
      alg = algo
    )
    alpha_old <- alpha
    alpha <- mle_out$solution
    
    alpha_mat <- matrix(
      alpha[(i - 1)*7 + wd],
      nrow = len_i,
      byrow = TRUE
    )
    
    diff_1 <- max(alpha - alpha_old)
    diff_2 <- max(A - A_old)
    diff_3 <- max(omega - omega_old)
    
    if (diff_1 < 1.0e-8 & diff_2 < 1.0e-8 & diff_3 < 1.0e-8) {
      break
    }
  }
  step2_converge <- run < max_iter
  
  # Forecasting next 2 weeks
  alpha_wd <- matrix(
    alpha, 
    nrow = len_i, 
    byrow = TRUE
  )
  wD <- wd_2[len_d]
  alpha_D <- alpha_mat[, len_d]
  u_D <- u_mat[, len_d]
  horizon_days <- horizon / len_t
  wd_h <- ((wD - 1 + 1:horizon_days) %% 7) + 1
  alpha_h <- alpha_wd[, wd_h]
  df_pred <- vector("list", horizon_days)
  
  for (k in 1:horizon_days) {
    ud_h <- alpha_h[, k] + expm::`%^%`(A, k) %*% (u_D - alpha_D)
    ud_h <- as.vector(ud_h)
    indx_f <- (1:len_i - 1) * 7 * len_t + (wd_h[k] - 1) * len_t
    indx_f <- map(indx_f, ~.x + 1:len_t) %>% unlist() 
    x_h <- rep(ud_h, each = len_t) * f_vec[indx_f]
    pred <- x_h^2 - 1/4
    
    df_pred[[k]] <- tibble(
      stream = rep(stream, each = len_t),
      h = rep(len_t * (k - 1) + 1:len_t, len_i),
      pred = pred
    )
  }
  df_pred <- bind_rows(df_pred)
  
  out <- list(
    df_pred = df_pred,
    step1_converge = step1_converge,
    step2_converge = step2_converge,
    params = list(u_vec = u_vec,
                  sigma = sigma,
                  f_vec = f_vec,
                  alpha = alpha,
                  A = A,
                  omega = omega)
  )
  return(out)
}



