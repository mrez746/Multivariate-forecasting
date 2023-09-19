initialize_fhat <- function(x, wd, t, hours,
                            weekdays = 1:7) {
  # Calculates f_hat from Eq (7) in Ye et al. (2019)
  # Step 1-a of The Iterative GLS Algorithm
  # Input:
  #       x (dbl vec) transformed call volumes 
  #       wd (int vec) weekday indices  
  #       t (int vec) intraday time interval indices  
  #       hours (int vec) unique indices of intraday time intervals, ordered  
  #       weekday (int vec) unique weekday indices, ordered  
  # Output:
  #       (dbl vec) of length |weekdays| * |hours|
  #           out[1]: f_hat wday[1] and hours[1], 
  #           out[2]: f_hat wday[1] and hours[2], etc.
  out <- rep(0, length(weekdays) * length(hours))
  
  for (i in 1:length(weekdays)) {
    
    for (j in 1:length(hours)) {
      indx_1 <- wd == weekdays[i]
      indx_2 <- t == hours[j]
      f_hat <- sum(x[indx_1 & indx_2]) / sum(x[indx_1])
      out[length(hours) * (i - 1) + j] <- f_hat
    }
  }
  return(out)
}

initialize_ud <- function(x, fhat, wd, d, t, hours,
                          weekdays = 1:7) {
  # Calculate u_d
  # Step 1-b of The Iterative GLS Algorithm
  # Inputs:
  #       x (dbl vec) transformed call volumes 
  #       fhat (dbl vec) estimated values of the f vector for 1 region/stream
  #         order: wday, time interval
  #       wd (int vec) weekday indices  
  #       t (int vec) intraday time interval indices  
  #       hours (int vec) unique indices of intraday time intervals, ordered  
  #       weekday (int vec) unique weekday indices, ordered  
  # Outputs:
  #       (tibble) w/ columns:
  #         $d
  #         $t
  #         $x
  #         $f
  #         $u_d
  #         $res
  indx <- length(hours) * (wd - 1) + t
  
  df_lm <- tibble(
    x = x, 
    f = fhat[indx],
    d = d,
    t = t
  )
  
  out <- df_lm %>% 
    group_by(d) %>% 
    nest() %>% 
    mutate(
      mod = map(data, ~lm(x ~ -1 + f, data = .x)),
      u_d = map_dbl(mod, ~coefficients(.x)["f"]),
      res = map(mod, ~.x$residuals)
    ) %>% 
    select(-mod) %>% 
    unnest(cols = c(data, res)) %>% 
    ungroup() %>% 
    arrange(d, t)
  return(out)
}

initialize_alpha <- function(len_i, i, wd, ud) {
  # Initializing for the procedure in Sec. 3.2.2 of Ye et al. (2019)
  # Inputs:
  #       len_i (int) number of time series
  #       i (int vec) specifying the time series of daily call volumes
  #       wd (list) specifying the weekday of call volumes 
  #       ud (list) the u_d values of daily call volumes
  # Outputs:
  #       (list) w/ components:
  #         $alpha (dbl vec) order: region ID, wday
  #         $sd_alpha (dbl vec) order: region ID, wday
  #         $alpha_mat (dbl mat) row # = region ID, col # = wday
  #         $u_mat (dbl mat) row # = region ID
  #         $n1_vec
  #         $n2_vec
  #         $one_list
  #         $wd
  
  alpha <- map2_dbl(
    rep(1:len_i, each = 7), 
    rep(1:7, len_i), 
    ~mean(ud[i == .x & wd == .y])
  )
  
  sd_alpha <- map2_dbl(
    rep(1:len_i, each = 7), 
    rep(1:7, len_i), 
    ~sd(ud[i == .x & wd == .y])
  )
  
  u_mat <- matrix(
    ud,
    nrow = len_i,
    byrow = TRUE
  )
  
  alpha_mat <- matrix(
    alpha[(i - 1)*7 + wd],
    nrow = len_i,
    byrow = TRUE
  )
  wd_2 <- wd[i == 1]
  n1_vec <- map_dbl(1:7, ~sum(wd_2[-length(wd_2)] == .x))
  n2_vec <- map_dbl(1:7, ~sum(wd_2[-1] == .x))
  one_list <- map(1:7, ~rep(1, n2_vec[.x]))
  
  out <- list(
    alpha = alpha,
    sd_alpha = sd_alpha,
    alpha_mat = alpha_mat,
    u_mat = u_mat,
    n1_vec = n1_vec,
    n2_vec = n2_vec,
    one_list = one_list
  )
  return(out)
}