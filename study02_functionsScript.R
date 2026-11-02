#' Pseudo weight estimation with Full Likelihood 
#' 
#' @param rep_id Number to identify Monte Carlo simulation
#' @param final_sample Dataframe of the sample
#' @param max_iter Maximum number of iterations for profiling the parameters
#' @param agb_wrapped Packed AGB raster object
#' @param covars_wrapped Packed Covariates raster object
#' @param aux_df Dataframe of the auxiliary data
#' @param true_mean Ture population mean
#' @return rep_id, naive mean of the sample, weighted mean of the sample, variance of the weighted mean, 
#'         lower and upper 95% confidence intervals

pesudo_weight_profileLik_pop <- function(rep_id, final_sample, max_iter, agb_wrapped, covars_wrapped, aux_df, ture_mean) {
  library(FNN)
  library(terra)
  library(alabama)
  library(Matrix)
  library(numDeriv)
  library(lpSolve)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)
  
  ## Extracting covariates at the non-probability sample locations
  covar_vals <- extract(covars, final_sample[, c("x", "y")])
  final_sample <- cbind(final_sample, covar_vals[, -1])
  final_sample_raw <- final_sample
  
  final_sample[colnames(final_sample)] <- scale(final_sample[colnames(final_sample)])
  
  cols <- c("x", "y", "NDVI", "EVI", "DEM", "SLOPE")
  
  # Z_matrix <- as.matrix(cbind(1, final_sample$x, final_sample$y, final_sample$NDVI, final_sample$EVI,
  #                             final_sample$DEM, final_sample$SLOPE, final_sample$AGB))
  
  Z_matrix <- cbind(1, as.matrix(final_sample[, cols]), AGB = final_sample$AGB)

  sample_size <- nrow(final_sample)
  
  mean_pop_NDVI <- mean(aux_df$NDVI)
  mean_pop_EVI <- mean(aux_df$EVI)
  mean_pop_DEM <- mean(aux_df$DEM)
  mean_pop_SLOPE <- mean(aux_df$SLOPE)
  mean_pop_x <- mean(aux_df$x)
  mean_pop_y <- mean(aux_df$y)
  
  
  # cols <- colnames(final_sample)[-3]
  mean_cols <- colMeans(final_sample_raw[, cols])
  std_cols <- apply(final_sample_raw[, cols], 2, sd)
  
  wt <- 1
  
  
  # Collecting the aux data for non sampled units
  sample_key <- paste(final_sample_raw$x, final_sample_raw$y)
  non_sample_key <- paste(aux_df$x, aux_df$y)
  
  non_sampled_units_raw <- aux_df[!(non_sample_key %in% sample_key), ]
  non_sampled_units <- non_sampled_units_raw
  
  sub_sample <- cbind(1, scale(non_sampled_units[, cols], center = mean_cols, scale = std_cols))
  sub_sample <- matrix(as.numeric(sub_sample), nrow = nrow(non_sampled_units))
  
  
  y_sample <- final_sample$AGB
  G_grid <- 4000
  
  Pi_Q <- function(gamma, Z_matrix, sub_sample, y_sample) {
    k <- ncol(sub_sample)
    list(
      Pi = as.vector(plogis(Z_matrix %*% gamma)),
      Q = plogis(-outer(as.vector(sub_sample %*% gamma[1:k]), gamma[k + 1] * y_sample, "+"))
    )
  }
  
  ## Likelihood for Gamma optimization, where p_i is assumed to be known
  neg_logLik_gamma <- function(gamma, Z_matrix, sub_sample, p, y_sample, sample_size, wt) {
    pq <- Pi_Q(gamma, Z_matrix, sub_sample, y_sample)
    
    Pi <- pq$Pi
    Q <- pq$Q
    
    if (any(Pi <= 0 | Pi >= 1)) return(1e10)
    
    h <- as.vector(Q %*% p)
    
    term1 <- sum(log(Pi))
    
    -(term1 + wt*sum(log(h)))
    
  }
  
  ## Gradient of the Gamma likelihood
  gradient_neg_logLik_gamma <- function(gamma, Z_matrix, sub_sample, p, y_sample, sample_size, wt) {
    pq <- Pi_Q(gamma, Z_matrix, sub_sample, y_sample)
    
    Pi <- pq$Pi
    Q <- pq$Q
    W <- Q * (1 - Q)
    
    h <- as.vector(Q %*% p)
    h1 <- -as.vector(W %*% p)
    hy <- -as.vector(W %*% (p*y_sample))
    
    term1 <- as.vector(crossprod(Z_matrix, 1 - Pi))
    term3 <- wt * c(as.vector(crossprod(sub_sample, h1/h)), sum(hy/h))
    
    -(term1 + term3)
  }
  
  ## Likelihood for p_i optimization, where Gamma is assumed to be known
  neg_logLik_p <- function(p, Q, wt) {
    
    if (any(p <= 0)) return(1e10)
    
    h <- as.vector(Q %*% p)
    term1 <- sum(log(p))
    
    -(term1 + wt*sum(log(h)))
  }
  
  ## Gradient of the likelihood
  gradient_neg_logLik_p <- function(p, Q, wt) {
    
    h <- as.vector(Q %*% p)
    
    -(1/p + wt*as.vector(crossprod(Q, 1/h)))
  }
  
  ## Equality constraints
  heq <- function(p, ...) {
    c(
      sum(p) - 1,
      sum(p * final_sample_raw$NDVI) - mean_pop_NDVI,
      sum(p * final_sample_raw$EVI) - mean_pop_EVI,
      sum(p * final_sample_raw$DEM) - mean_pop_DEM,
      sum(p * final_sample_raw$SLOPE) - mean_pop_SLOPE,
      sum(p * final_sample_raw$x) -  mean_pop_x,
      sum(p * final_sample_raw$y) - mean_pop_y
    )
  }
  
  ## Jacobian matrix for equality constraints
  heq.jac <- function(p, ...) {
    rbind(
      rep(1, length(p)),
      final_sample_raw$NDVI,
      final_sample_raw$EVI,
      final_sample_raw$DEM,
      final_sample_raw$SLOPE,
      final_sample_raw$x,
      final_sample_raw$y
    )
  }
  
  make_constraints <- function(p) {
    Cp <- heq.jac(p)
    b <- as.vector(Cp %*% p) - heq(p)
    m <- rowMeans(Cp[-1, , drop = FALSE])
    s <- apply(Cp[-1, , drop = FALSE], 1, sd)
    
    list(C = rbind(Cp[1, ], (Cp[-1, , drop = FALSE] - m)/s), 
         b = c(b[1], (b[-1] - m)/s))
    
  }
  
  ## Inequality constraints
  hin <- function(p, ...) {
    p - 1e-8
  }
  
  ## Jacobian matrix for inequality constraints
  hin.jac <- function(p, ...) diag(length(p))
  
  
  ## Full likelihood
  full_logLik <- function(gamma, p) {
    pq <- Pi_Q(gamma, Z_matrix, sub_sample, y_sample)
    h <- as.vector(pq$Q %*% p)
    
    (sum(log(pq$Pi)) + sum(log(p)) + sum(log(h)))
  }
  
  
  ## Newton Method Function
  newton_method <- function(p, Q, C, b, wt, max_iter = 20, tol = 1e-08){
    m <- nrow(C)
    
    r <- as.vector(C %*% p) - b
    p <- p - as.vector(crossprod(C, solve(tcrossprod(C), r)))
    stopifnot(all(p > 0))
    F0 <- neg_logLik_p(p, Q, wt)
    
    for (i in 1:max_iter) {
      h <- as.vector(Q %*% p)
      g <- -(1/p + wt*as.vector(crossprod(Q, 1/h)))
      H <- diag(1/p^2, length(p)) + wt * crossprod(Q/h)
      sol <- solve(rbind(cbind(H, t(C)), cbind(C, matrix(0, m, m))), c(-g, rep(0, m)))
      dp <- sol[seq_along(p)]
      dec <- -sum(g * dp)
      
      t <- 1
      if (any(dp < 0)) t <- min(1, 0.95 * min(-p[dp < 0]/dp[dp < 0]))
      repeat {
        pn <- p + t * dp
        Fn <- neg_logLik_p(pn, Q, wt)
        if (is.finite(Fn) && Fn <= F0 - 1e-04*t*dec) break
        t <- t /2
        if (t < 1e-12) break
      }
      p <- pn
      F0 <- Fn
      if (dec/2 < tol) break
    }
    list(p = p, value = F0, iter = i, decrement = dec)
  }
  
  ## Initialize gamma and p values for optimization
  gamma_init <- rnorm(ncol(Z_matrix), mean = 0, sd = 0.01)
  p_init <- rep(1/sample_size, sample_size)
  
  threshold <- 1e-06      ## Threshold for stopping when convergence is achieved
  max_iter <- max_iter    ## Max iterations for the Profile Lik to run
  
  ## Intial p's that satisfy all the constraints, linear programming
  interior_start <- function(C, b) {
    n <- ncol(C)
    m <- nrow(C)
    obj  <- c(rep(0, n), 1)
    A_eq <- cbind(C, 0)
    A_ge <- cbind(diag(n), -1)
    sol  <- lp("max", obj, rbind(A_eq, A_ge),
               c(rep("=", m), rep(">=", n)), c(b, rep(0, n)))
    list(p = sol$solution[1:n], t = sol$solution[n + 1], status = sol$status)
  }
  
  con <- make_constraints(p_init)
  st  <- interior_start(con$C, con$b)
  
  p_init <- st$p
  
  tol_mu <- 1e-3
  tol_LL <- 1e-3
  tol_gamma <- 1e-4
  tol_p <- 1e-5
  
  mu_hat <- numeric(max_iter)
  LL_full <- numeric(max_iter)
  
  ## Loop for Profile Likelihood
  for (iter in 1:max_iter) {
    
    ## Gamma optimization
    gamma_fit <- optim(
      par = gamma_init,
      fn = neg_logLik_gamma,
      gr = gradient_neg_logLik_gamma,
      Z_matrix = Z_matrix,
      sub_sample = sub_sample,
      y_sample = y_sample,
      p = p_init,
      sample_size = sample_size,
      wt = wt,
      method = "BFGS"
    )
    
    ## Collecting the new gamma parameter values
    gamma_new <- gamma_fit$par
    
    ## p Optimization
    pq <- Pi_Q(gamma_new, Z_matrix, sub_sample, y_sample)
    p_fit <- newton_method(p_init, pq$Q, con$C, con$b, wt)
    
    ## Collecting the new p parameter values
    p_new <- p_fit$p
    
    if (any(!is.finite(p_new))) {
      return(c(rep = rep_id, naive_mean = mean(final_sample_raw$AGB), mu_y_EL = NA))
    }
    
    p_new <- pmax(p_new, 1e-8)

    mu_hat[iter] <- sum(p_new * final_sample_raw$AGB)
    LL_full[iter] <- full_logLik(gamma_new, p_new)
    
    ## Checking if the change compared to previous iteration is less than set threshold
    if(iter > 1) {
      delta_mu <- abs(mu_hat[iter] - mu_hat[iter -1])
      delta_LL <- abs(LL_full[iter] - LL_full[iter - 1])
      delta_p <- max(abs(p_new - p_init))
      delta_gamma <- max(abs(gamma_new - gamma_init))
      if (delta_mu < tol_mu && delta_LL < tol_LL && delta_gamma < tol_gamma && delta_p < tol_p) {
        break
      }
    }
    
    ## Re-initalizing gamma and p values, as the new parameter values, for the next iteration
    gamma_init <- gamma_new
    p_init <- p_new
    
  }
  
  ## Naive mean of the sample
  naive_mean <- mean(final_sample_raw$AGB)
  
  ## Mean estimation using the weights (p_i)
  mu_y_EL <- sum(p_new * final_sample_raw$AGB)
  
  ## Collecting the optimal gamma and p_i
  optim_params <- c(gamma_new, p_new)
  n_g <- length(gamma_new)
  n_p <- length(p_new)
  
  
  ## Variance computation - Taylor series apporixmation, Inverse of 2nd Derivative of the LL, evaluated at optimal parameters
  k <- ncol(sub_sample)
  Q <- plogis(-outer(as.vector(sub_sample %*% gamma_new[1:k]), gamma_new[k + 1] * y_sample, "+"))
  h <- as.vector(Q %*% p_new)
  
  Hpp <- (diag(1/p_new^2) + crossprod(Q/h))
  
  Hgg <- jacobian(function(g) gradient_neg_logLik_gamma(g, Z_matrix, sub_sample, p_new, y_sample, sample_size, wt), 
                  x = gamma_new, method = "simple")
  
  # Hgp <- jacobian(function(p) gradient_neg_logLik_gamma(gamma_new, Z_matrix, sub_sample, p, y_sample, sample_size, wt), 
  #                             x = p_new, method = "simple")
  # 
  # Hpg <- t(Hgp)
  
  Hpg <- jacobian(function(g) gradient_neg_logLik_p(p_new, Pi_Q(g, Z_matrix, sub_sample, y_sample)$Q, wt),
                  x = gamma_new, method = "simple")
  
  Hgp <- t(Hpg)
  
  H_block <- rbind(cbind(Hgg, Hgp), cbind(Hpg, Hpp))
  
  H_block <- (H_block + t(H_block))/2
  
  Cp <- heq.jac(p_new)
  Cs <- rbind(Cp[1, ], t(scale(t(Cp[-1, ]))))
  
  C <- cbind(matrix(0, nrow(Cs), n_g), Cs)
  
  # C <- cbind(matrix(0, nrow(heq.jac(p_new)), n_g), heq.jac(p_new))
  
  K <- rbind(cbind(H_block, t(C)),
             cbind(C, matrix(0, nrow(C), nrow(C))))
  
  V <- solve(K)[1:(n_g + n_p), 1: (n_g + n_p)]
  
  g_vec <- c(rep(0, n_g), final_sample_raw$AGB)
  var_mu <- drop(t(g_vec) %*% V %*% g_vec)
  
  ci_u <- mu_y_EL + 1.96 * sqrt(var_mu)
  ci_l <- mu_y_EL - 1.96 * sqrt(var_mu)
  
  return(c(rep = rep_id, naive_mean = naive_mean, mu_y_EL = mu_y_EL, var_mu = var_mu, ci_l = ci_l, ci_u = ci_u))
}
