library(sp)
library(sf)
library(ggplot2)
library(tidyr)
library(dplyr)
library(terra)
library(spcosa)
library(FNN)
library(future.apply)
library(spatstat.geom)
library(spatstat.random)
library(alabama)
library(UPSvarApprox)
library(samplingVarEst)
library(moments)
library(gstat)
library(mgcv)
library(clhs)
library(Matrix)
library(numDeriv)
library(lpSolve)


## Load the Target Variable
# agb <- rast("../study02_data/study02_agb_simulated_exp_NDVI_EVI_expCorrResid150.tif")
# agb <- rast("../study02_results/study02_agg/study02_agb_simulated_agg_exp_NDVI_EVI_expCorrResid150.tif")
agb <- rast("../study02_results/study02_agg/study02_agb_simulated_agg_exp_NDVI_EVI_AGB_expCorrResid180.tif")
# covars <- rast("../study02_data/study02_env_covariates#4.tif")
covars <- rast("../study02_data/study02_env_covariates#4_filtered.tif")
plot(agb)


## Prepare Spatial DF of the target variable
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
coordinates(agb_df) <- ~x+y
gridded(agb_df) <- TRUE
agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)

## Collecting Covariate names
covars_names <- names(covars)

# agb_cci <- rast("../study02_data/study02_agb_CCI_filtered.tif")
# agb_cci_df <- as.data.frame(agb_cci, xy = TRUE, na.rm = TRUE)
valid_cells <- ifel(is.na(sum(covars)), NA, 1)

covars_agg <- aggregate(covars * valid_cells, fact = 10, fun = "mean", na.rm = TRUE)
# agb_agg <- aggregate(agb_cci * valid_cells, fact = 10, fun = "mean", na.rm = TRUE)

names(covars_agg) <- names(covars)
# names(agb_agg) <- "AGB"

# agb_cci <- agb_agg
covars <- covars_agg

covars <- c(covars, agb_cci)
aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)

## Standardize AGB values
m <- global(agb, "mean", na.rm = TRUE)[1, 1]
s <- global(agb, "sd", na.rm = TRUE) [1, 1]
agb_std <- (agb - m)/s

covars_m <- global(covars, "mean", na.rm = TRUE)$mean
covars_sd <- global(covars, "sd", na.rm = TRUE)$sd
covars_std <- (covars - covars_m)/covars_sd

## Wrap the Spatraster - for parallel processing
covars_wrapped <- wrap(covars)
agb_wrapped <- wrap(agb)
agb_std_wrapped <- wrap(agb_std)
covars_std_wrapped <- wrap(covars_std)

## Function for generating Preferential Samples - Strength is controlled by var - Beta
pref_sample <- function(alpha = 0, beta = 1, sample_size, agb_wrapped, agb_std_wrapped) {
  
  agb <- unwrap(agb_wrapped)
  agb_std <- unwrap(agb_std_wrapped)
  
  ## Sampling Intensity - based on the Response
  lambda <- exp(alpha + beta * agb_std)
  v <- values(lambda)
  expected_n <- sample_size
  
  scale_factor <- expected_n/sum(v, na.rm = TRUE)
  lambda_scaled <- lambda * scale_factor
  lam_values <- values(lambda_scaled)
  counts <- rep(NA, length(lam_values))
  valid <- !is.na(lam_values)
  
  counts[valid] <- rpois(sum(valid), lam_values[valid])
  
  cells_positive <- rep(which(!is.na(counts) & counts > 0))
  
  cells <- unlist(
    lapply(cells_positive, function(i) rep(i, counts[i]))
  )
  
  
  ## Extracting XY and AGB values for the points
  xy <- xyFromCell(agb, cells)
  
  final_sample <- data.frame(
    x = xy[, 1],
    y = xy[, 2]
  )
  
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}

pref_resp_covar_sample <- function(alpha = 0, beta0 = 1, beta1 = 1, sample_size, agb_wrapped, agb_std_wrapped, covars_std_wrapped) {
  
  agb <- unwrap(agb_wrapped)
  agb_std <- unwrap(agb_std_wrapped)
  covars_std <- unwrap(covars_std_wrapped)
  
  
  ## Sampling Intensity - based on the Response
  lambda <- exp(alpha + beta0 * agb_std - beta1 * covars_std[[3]] - beta1 * covars_std[[4]])
  v <- values(lambda)
  expected_n <- sample_size
  
  scale_factor <- expected_n/sum(v, na.rm = TRUE)
  lambda_scaled <- lambda * scale_factor
  lam_values <- values(lambda_scaled)
  counts <- rep(NA, length(lam_values))
  valid <- !is.na(lam_values)
  
  counts[valid] <- rpois(sum(valid), lam_values[valid])
  
  cells_positive <- rep(which(!is.na(counts) & counts > 0))
  
  cells <- unlist(
    lapply(cells_positive, function(i) rep(i, counts[i]))
  )
  
  
  ## Extracting XY and AGB values for the points
  xy <- xyFromCell(agb, cells)
  
  final_sample <- data.frame(
    x = xy[, 1],
    y = xy[, 2]
  )
  
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}


## Function for generating strongly clustered samples - strength depends on param - s
strongly_clustered_sample <- function(s = 0.9, sample_size, agb_wrapped, strata_IDs) {
  
  agb <- unwrap(agb_wrapped)
  
  ## Randomly drawing 2 stratas - one with 10 stratas, and one with 90 stratas
  cluster_strata <- sample(strata_IDs, 10)
  other_strata <- setdiff(strata_IDs, cluster_strata)
  
  
  ## Drawing 90% sample size from the small group of stratas, and 10% from the large group of stratas
  sample_size <- 500
  n_cluster <- round(s*sample_size)
  n_other <- sample_size - n_cluster
  
  
  ## Extracting all the pixels from the both strata
  cluster_df <- strata_map_df[strata_map_df$stratumId %in% cluster_strata, ]
  other_df <- strata_map_df[strata_map_df$stratumId %in% other_strata, ]
  
  
  ## Drawing samples from those two strata
  cluster_sample <- cluster_df[sample(1:nrow(cluster_df), n_cluster), ]
  other_sample <- other_df[sample(1:nrow(other_df), n_other), ]
  
  
  ## Extracting target variable at these locations
  final_sample <- rbind(data.frame(x = cluster_sample$x, y = cluster_sample$y, ID = "cluster_stratum"),
                        data.frame(x = other_sample$x, y = other_sample$y, ID = "other_stratum"))
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  final_sample <- final_sample[, c(1, 2, 4)]
  
  return(final_sample)
}

## Strongly Clustered Samples done with Compact GeoStrata DOI: 10.1016/j.ecoinf.2022.101665
set.seed(123)
n_strata <- 100
stratification <- stratify(agb_df, nStrata = n_strata, nTry = 1)
strata_map <- as(stratification, "SpatialPixelsDataFrame")
strata_map_df <- as.data.frame(strata_map)
strata_IDs <- unique(strata_map_df$stratumId)


## Function for generating Conditioned Latin Hypercube Sampling
clhs_sample <- function(sample_size = 500, agb_wrapped, covars_df = covars_df_new){
  
  agb <- unwrap(agb_wrapped)
  
  clhs_result <- clhs(covars_df, size = sample_size, cost = "DistRoads", iter = 10000, simple = FALSE, progress = TRUE)
  
  sample_data <- clhs_result$sampled_data
  
  final_sample <- data.frame(sample_data[, c("x", "y")], row.names = NULL)
  
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
  
}

accessibility_sample <- function(sample_size = 500, agb_wrapped, covars_df = covars_df_new, a1 = 1, a2 = 0.5, a3 = 0.5) {
  
  agb <- unwrap(agb_wrapped)
  access_df <- covars_df
  access_df$DistRoads_s <- scale(access_df$DistRoads)
  access_df$DEM_s <- scale(access_df$DEM)
  access_df$SLOPE_s <- scale(access_df$SLOPE)
  
  access_df$Cost <- a1*access_df$DistRoads_s + a2*access_df$DEM_s + a3*access_df$SLOPE_s
  
  prob <- exp(-access_df$Cost)
  prob <- prob/sum(prob)
  
  id <- sample(seq(nrow(access_df)), size = sample_size, replace = FALSE, prob = prob)
  
  sample_data <- access_df[id, ]
  
  final_sample <- data.frame(
    x = sample_data$x,
    y = sample_data$y,
    row.names = NULL
  )
  
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}

roads <- st_read("../study02_data/study02_roadsClipped.shp")
roads <- st_transform(roads, crs(covars))
roads_vect <- vect(roads)

dist_to_roads <- distance(covars[[1]], roads_vect)

covars_stack <- c(covars, dist_to_roads)
covars_df_new <- as.data.frame(covars_stack, xy = TRUE, na.rm = TRUE)
colnames(covars_df_new)[7] <- "DistRoads"

final_sample <- accessibility_sample(500, agb_wrapped, covars_df_new, 5, 0, 0)
summary(final_sample)
## Empirical Likelihood implementation, with Profile Likelihood
pesudo_weight_profileLik <- function(rep_id, final_sample, max_iter, agb_wrapped, covars_wrapped, aux_df) {
  library(FNN)
  library(terra)
  library(alabama)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)
  
  ## Extracting covariates at the non-probability sample locations
  covar_vals <- terra::extract(covars, final_sample[, c("x", "y")])
  final_sample <- cbind(final_sample, covar_vals[, -1])
  final_sample_raw <- final_sample

  final_sample[colnames(final_sample)] <- scale(final_sample[colnames(final_sample)])
  
  n <- nrow(final_sample)
  
  ## Mean population estimates for auxiliary data
  mean_pop_NDVI <- mean(aux_df$NDVI)
  mean_pop_EVI <- mean(aux_df$EVI)
  mean_pop_DEM <- mean(aux_df$DEM)
  mean_pop_SLOPE <- mean(aux_df$SLOPE)
  
  
  ## Matrix of covariates for estimating inclusion probabilities
  Z_matrix <- as.matrix(cbind(1, final_sample$x, final_sample$y, final_sample$NDVI, final_sample$EVI, 
    final_sample$DEM, final_sample$SLOPE, final_sample$AGB))
  
  
  ## Inclusion probability function, logistic regression cdf, gamma - is the nuisance parameter
  P_i <- function(gamma, Z_matrix) {
    plogis(Z_matrix %*% gamma)
  }
  
  ## Likelihood for Gamma optimization, where p_i is assumed to be known
  neg_logLik_gamma <- function(gamma, Z_matrix, p, n){
    
    if (any(p <= 0)) return(1e10)
    
    Pi <- P_i(gamma, Z_matrix)
    
    if (any(Pi <= 0 | Pi >= 1)) return(1e10)
    
    term1 <- sum(log(Pi))
    term2 <- n * log(sum(Pi * p))
    
    -(term1 - term2)
  }
  
  ## Likelihood for p_i optimization, where Gamma is assumed to be known
  neg_logLik_p <- function(p, n, gamma, Z_matrix) {
    
    if (any(p <= 0)) return(1e10)
    
    Pi <- P_i(gamma, Z_matrix)
    
    term1 <- n * log(sum(p*Pi))
    term2 <- sum(log(p))
    -(term2 - term1)
  }
  
  ## Gradient of the likelihood
  gradient_neg_logLik_p <- function(p, n, gamma, Z_matrix) {
    
    Pi <- P_i(gamma, Z_matrix)
    
    S <- sum(p*Pi)
    
    -(1/p - n*Pi/S)
  }
  
  ## Equality constraints
  heq <- function(p, ...) {
    c(
      sum(p) - 1,
      sum(p * final_sample_raw$NDVI) - mean_pop_NDVI,
      sum(p * final_sample_raw$EVI) - mean_pop_EVI,
      sum(p * final_sample_raw$DEM) - mean_pop_DEM,
      sum(p * final_sample_raw$SLOPE) - mean_pop_SLOPE
    )
  }
  
  ## Jacobian matrix for equality constraints
  heq.jac <- function(p, ...) {
    rbind(
      rep(1, length(p)),
      final_sample_raw$NDVI,
      final_sample_raw$EVI,
      final_sample_raw$DEM,
      final_sample_raw$SLOPE
    )
  }
  
  ## Inequality constraints
  hin <- function(p, ...) {
    p - 1e-8
  }
  
  ## Jacobian matrix for inequality constraints
  hin.jac <- function(p, ...) diag(length(p))
  
  ## Initialize gamma and p values for optimization
  gamma_init <- rnorm(ncol(Z_matrix), mean = 0, sd = 0.01)
  p_init <- rep(1/n, n)

  threshold <- 1e-06      ## Threshold for stopping when convergence is achieved
  max_iter <- max_iter    ## Max iterations for the Profile Lik to run

  ## Loop for Profile Likelihood
  for (iter in 1:max_iter) {

    ## Gamma optimization
    gamma_fit <- optim(
      par = gamma_init,
      fn = neg_logLik_gamma,
      Z_matrix = Z_matrix,
      p = p_init,
      n = n,
      method = "BFGS"
    )
    
    ## Collecting the new gamma parameter values
    gamma_new <- gamma_fit$par

    ## p Optimization
    p_fit <- auglag(
      par = p_init,
      fn = neg_logLik_p,
      gr = gradient_neg_logLik_p,
      gamma = gamma_new,
      Z_matrix = Z_matrix,
      n = n,
      heq = heq,
      heq.jac = heq.jac,
      hin = hin,
      hin.jac = hin.jac,
      control.outer = list(trace = FALSE)
    )

    ## Collecting the new p parameter values
    p_new <- p_fit$par

    ## Checking if the change compared to previous iteration is less than set threshold
    delta_p <- max(abs(p_new - p_init))
    delta_gamma <- max(abs(gamma_new - gamma_init))

    ## If threshold reached, break the loop
    if (delta_p < threshold && delta_gamma < threshold) {
      break
    }

    ## Re-initalizing gamma and p values, as the new parameter values, for the next iteration
    gamma_init <- gamma_new
    p_init <- p_new

  }

  ## Naive mean of the sample
  naive_mean <- mean(final_sample_raw$AGB)

  ## Mean estimation using the weights (p_i)
  mu_y_EL <- sum(p_new * final_sample_raw$AGB)
  
  return(c(rep = rep_id, naive_mean = naive_mean, mu_y_EL = mu_y_EL))
}

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

    mu_PL[iter] <- sum(p_new * final_sample_pref_raw$AGB)
    LL_full[iter] <- full_logLik(gamma_new, p_new)
    
    ## Checking if the change compared to previous iteration is less than set threshold
    if(iter > 1) {
      delta_mu <- abs(mu_hat[iter] - mu_PL[iter -1])
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

pseudo_weight_knn_EL_clac <- function(rep_id, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean) {
  
  library(FNN)
  library(terra)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)
  
  extracted <- terra::extract(covars, final_sample[, c("x", "y")])
  final_sample <- cbind(final_sample, extracted[, covars_names])
  
  raw_covars <- final_sample[, covars_names]
  colnames(raw_covars) <- paste0(covars_names, "_raw")
  
  
  ## Normalizing Covariates for KNN
  norm_func <- function(df, colnames, min_vals = NULL, max_vals = NULL){
    for (i in colnames){
      df[[i]] <- (df[[i]] - min_vals[i])/(max_vals[i] - min_vals[i])
    }
    return(df)
  }
  
  B_min <- sapply(final_sample[covars_names], min, na.rm = TRUE)
  B_max <- sapply(final_sample[covars_names], max, na.rm = TRUE)
  
  
  aux_norm_df <- norm_func(aux_df, covars_names, B_min, B_max)
  final_sample_norm_df <- norm_func(final_sample, covars_names, B_min, B_max)
  final_sample_norm_df <- cbind(final_sample_norm_df, raw_covars)
  
  
  ## Creating a membership indicator for units in A
  aux_norm_df$mem_indicator <- paste(aux_norm_df$x, aux_norm_df$y) %in% paste(final_sample_norm_df$x, final_sample_norm_df$y)
  match_row <- match(paste(aux_norm_df$x, aux_norm_df$y), paste(final_sample_norm_df$x, final_sample_norm_df$y))
  
  
  ## Matrices for KNN
  query_matrix <- as.matrix(aux_norm_df[, covars_names])
  data_matrix <- as.matrix(final_sample_norm_df[, covars_names])
  
  
  ## KNN
  k <- 2
  knn_results <- get.knnx(data_matrix, query_matrix, k + 1) ## computes the ids of the closest K neighbors, along with distances, for units in query matrix
  
  
  ## Getting AGB values using the closest K neighbors, for B units, dropping its true y, so that there is a prediction there
  agb_impute <- numeric(nrow(aux_norm_df))
  
  for (i in seq_len(nrow(aux_norm_df))) {
    id <- knn_results$nn.index[i, ]
    
    if (!is.na(match_row[i])) {
      id <- id[id != match_row[i]][1:k]
    } else {
      id <- id[1:k]
    }
    
    agb_impute[i] <- mean(final_sample_norm_df$AGB[id])
  }
  
  
  aux_norm_df$agb_impute <- agb_impute
  
  formula <- as.formula(paste("mem_indicator ~ x + y + agb_impute + ", paste(covars_names, collapse = " + ")))
  
  ## Logistic regression of the membership indicator
  sampling_model <- glm(formula, family = binomial(link = "logit"), data = aux_norm_df)
  # sampling_model <- glm(mem_indicator ~ NDVI + EVI + DEM + SLOPE + agb_impute, family = binomial(link = "logit"), data = aux_norm_df)
  
  non_prob_sample_pred <- final_sample_norm_df
  non_prob_sample_pred$agb_impute <- final_sample_norm_df$AGB
  
  
  ## Predicting inclusion probabilities & pseudo weights for the units of the non-probability sample
  final_sample_norm_df$inc_probs <- predict(sampling_model, newdata = non_prob_sample_pred, type = "response")
  
  final_sample_norm_df$pseudo_weights <- 1/final_sample_norm_df$inc_probs
  
  
  ## HT Estimator
  n <- nrow(aux_norm_df)
  
  naive_mean <- mean(final_sample_norm_df$AGB)
  pseudo_weight_mean <- (1/n) * sum(final_sample_norm_df$pseudo_weights * final_sample_norm_df$AGB)
  
  
  ## Empirical Likelihood
  sample_size <- nrow(final_sample_norm_df)
  pseudo_weights <- final_sample_norm_df$pseudo_weights
  
  mean_NDVI_pop <- mean(aux_df$NDVI)
  mean_EVI_pop <- mean(aux_df$EVI)
  mean_DEM_pop <- mean(aux_df$DEM)
  mean_SLOPE_pop <- mean(aux_df$SLOPE)
  
  ## Objective function of the EL
  neg_logLik <- function(p) {
    if (any(p <= 0)) return(1e10)
    term1 <- sample_size * log(sum(p*pseudo_weights))
    term2 <- sum(log(p))
    -(term1 + term2)
  }
  
  gr_logLik <- function(p) {
    S <- sum(p*pseudo_weights)
    -((sample_size*pseudo_weights)/S + 1/p)
  }
  
  
  ## Equality constraints
  heq <- function(p) {
    c(
      sum(p) - 1,
      sum(p * final_sample_norm_df$NDVI_raw) - mean_NDVI_pop,
      sum(p * final_sample_norm_df$EVI_raw) - mean_EVI_pop,
      sum(p * final_sample_norm_df$DEM_raw) - mean_DEM_pop,
      sum(p * final_sample_norm_df$SLOPE_raw) - mean_SLOPE_pop
    )
  }
  
  heq.jac <- function(p) {
    rbind(
      rep(1, sample_size),
      final_sample_norm_df$NDVI_raw,
      final_sample_norm_df$EVI_raw,
      final_sample_norm_df$DEM_raw,
      final_sample_norm_df$SLOPE_raw
    )
  }
  
  ## Inequality constraints
  hin <- function(p) {
    p - 1e-8
  }
  
  hin.jac <- function(p) diag(sample_size)
  
  ## Initializing point masses
  p_init <- rep(1/sample_size, sample_size)
  
  
  ## Augmented Lagrangian - for constrained optimization
  result_EL <- auglag(par = p_init, fn = neg_logLik, gr = gr_logLik, heq = heq, heq.jac = heq.jac,
                      hin = hin, hin.jac = hin.jac, control.outer = list(trace = FALSE))
  
  p_hat <- result_EL$par  #point masses
  
  ## HT estimator using the weights
  mu_y_EL <- sum(p_hat * final_sample_norm_df$AGB)
  
  
  # return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, mu_y_EL = mu_y_EL,
  #             final_sample_norm_df = final_sample_norm_df))
  return(c(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, mu_y_EL = mu_y_EL))
  
}


## Parallel computing: For repeating the estimation process over MC simulations
plan(multisession, workers = parallel::detectCores(logical = FALSE) - 1)

true_mean <- mean(agb_df_plot$AGB)
mc_reps <- 50                       ## Monte Carlo simulations
sample_size <- 500

set.seed(123)
results_list <- future_lapply(
  1:mc_reps,
  function(r) {
    final_sample <- pref_sample(alpha = 0, beta = 0.5, sample_size = sample_size, agb_wrapped, agb_std_wrapped)
    # final_sample <- pref_resp_covar_sample(alpha = 0, beta0 = 1, beta1 = 0.5, sample_size = 500, agb_wrapped, agb_std_wrapped, covars_std_wrapped)
    # final_sample <- strongly_clustered_sample(s = 0.9, sample_size, agb_wrapped, strata_IDs)
    # final_sample <- clhs_sample(sample_size = sample_size, agb_wrapped, covars_df = covars_df_new)
    # result <- pesudo_weight_profileLik(r, final_sample, max_iter = 10, agb_wrapped, covars_wrapped, aux_df)
    # result <- pseudo_weight_knn_EL_clac(r, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
    # result <- pesudo_weight_profileLik_pop(r, final_sample, max_iter = 5, agb_wrapped, covars_wrapped, aux_df, true_mean)
    result <- tryCatch(
      pesudo_weight_profileLik_pop(r, final_sample, max_iter = 15, agb_wrapped, covars_wrapped, aux_df, ture_mean),
      error = function(e) {
        cat(sprintf("MC iter %d failed: %s\n", r, conditionMessage(e)),
            file = "progress_log.txt", append = TRUE)
        c(rep = r, naive_mean = NA, mu_y_EL = NA)
      })

    cat(sprintf("MC iter %d done | naive_mean = %.3f | mu_y_EL = %.3f | var = %.3f | ci_l = %.3f | ci_u = %.3f\n",
                r, result["naive_mean"], result["mu_y_EL"], result["var_mu"], result["ci_l"], result["ci_u"]),
        file = "progress_log.txt", append = TRUE)
    result
    
    # cat(sprintf("MC iter %d done\n", r), file = "progress_log.txt", append = TRUE)
    # final_sample
    },
  future.seed = TRUE
)

plan(sequential)

results <- as.data.frame(do.call(rbind, results_list))

write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_prefSample05_NewLik_var.csv", row.names = FALSE)
results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample05_NewLik_var.csv")

head(results)

plot_df <- data.frame(
  estimator = rep(c("naive_mean", "mu_y_EL"), each = nrow(results)),
  estimate = c(results$naive_mean, results$mu_y_EL)
)


ggplot(plot_df, aes(x = estimator, y = estimate, fill = estimator)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("naive_mean" = "Naive", "mu_y_EL" = "PL Wght")) +
  labs(title = "Naive, EL Mean Estimator",
       subtitle = paste("Dashed line = true mean:", round(true_mean, 2)),
       y = "AGB estimate", x = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

results$coverage <- results$ci_l <= true_mean & results$ci_u >= true_mean

accuracy_metrics <- function(df, pop_mean) {
  
  bias_pl <- mean(df$mu_y_EL, na.rm = TRUE) - pop_mean
  sd_pl <- sd(df$mu_y_EL, na.rm = TRUE)
  rmse_pl <- sqrt(mean((df$mu_y_EL - pop_mean)^2, na.rm = TRUE))
  coverage_pl <- mean(df$coverage, na.rm = TRUE)* 100
  
  return(data.frame(
    "PL_Bias" = round(bias_pl, 2),
    "PL_Std" = round(sd_pl, 2),
    "PL_RMSE" = round(rmse_pl, 2),
    "PL_Coverage (%)" = round(coverage_pl, 2),
    check.names = FALSE
  ))
}

results_accuracy <- accuracy_metrics(results, true_mean)
results_accuracy
####################TRUE#####################################################################
## Profile Likelihood Testing
set.seed(12345)
final_sample_pref <- pref_sample(0, 1, 500, agb_wrapped, agb_std_wrapped)
# final_sample_pref <- pref_resp_covar_sample(alpha = 0, beta0 = 1, beta1 = 0.5, sample_size = 500, agb_wrapped, agb_std_wrapped, covars_std_wrapped)
# final_sample_pref <- strongly_clustered_sample(s = 0.9, 500, agb_wrapped, strata_IDs)
# final_sample_pref <- clhs_sample(500, agb_wrapped, covars_df_new)
# final_sample_pref <- results_list[[1]]

final_sample_pref_sf <- st_as_sf(final_sample_pref, coords = c("x", "y"), crs = crs(covars))

covar_vals <- extract(covars, final_sample_pref[, c("x", "y")])
final_sample_pref <- cbind(final_sample_pref, covar_vals[, -1])
final_sample_pref_raw <- final_sample_pref

final_sample_pref[colnames(final_sample_pref)] <- scale(final_sample_pref[colnames(final_sample_pref)])

# cols <- c("x", "y", "NDVI", "EVI", "DEM", "SLOPE")
cols <- c("x", "y", "NDVI", "EVI", "DEM", "SLOPE")

# Z_matrix <- as.matrix(cbind(1, final_sample_pref$x, final_sample_pref$y, final_sample_pref$NDVI, final_sample_pref$EVI,
#   final_sample_pref$DEM, final_sample_pref$SLOPE, final_sample_pref$AGB))

Z_matrix <- cbind(1, as.matrix(final_sample_pref[, cols]), AGB = final_sample_pref$AGB)


sample_size <- nrow(final_sample_pref)

mean_pop_NDVI <- mean(aux_df$NDVI)
mean_pop_EVI <- mean(aux_df$EVI)
mean_pop_DEM <- mean(aux_df$DEM)
mean_pop_SLOPE <- mean(aux_df$SLOPE)
mean_pop_x <- mean(aux_df$x)
mean_pop_y <- mean(aux_df$y)


mean_cols <- colMeans(final_sample_pref_raw[, cols])
std_cols <- apply(final_sample_pref_raw[, cols], 2, sd)


# 
# block_size <- 4000
# 
# sub_sample <- X_aux[sample(nrow(non_sampled_units), 5000), ]

# wn <- nrow(aux_norm_df)t <- nrow(non_sampled_units)/5000
wt <- 1


sample_key <- paste(final_sample_pref_raw$x, final_sample_pref_raw$y)
non_sample_key <- paste(aux_df$x, aux_df$y)

non_sampled_units_raw <- aux_df[!(non_sample_key %in% sample_key), ]
non_sampled_units <- non_sampled_units_raw

# sub <- sample(nrow(non_sampled_units), 5000)
# sub_sample <- cbind(1, scale(non_sampled_units[sub, cols], center = mean_cols, scale = std_cols))
# sub_sample <- matrix(as.numeric(sub_sample), nrow = 5000)


sub_sample <- cbind(1, scale(non_sampled_units[, cols], center = mean_cols, scale = std_cols))
sub_sample <- matrix(as.numeric(sub_sample), nrow = nrow(non_sampled_units))


y_sample <- final_sample_pref$AGB
G_grid <- 4000

## Inclusion probability function, logistic regression cdf, gamma - is the nuisance parameter
P_i <- function(gamma, Z_matrix) {
  plogis(Z_matrix %*% gamma)
}

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
  # term2 <- sample_size * log(sum(Pi * p))

  -(term1 + wt*sum(log(h)))

}

# neg_logLik_gamma <- function(gamma, Z_matrix, sub_sample, p, y_sample, sample_size, wt) {
#   Pi <- as.vector(plogis(Z_matrix %*% gamma))
# 
#   if (any(Pi <= 0 | Pi >= 1)) return(1e10)
# 
#   a <- as.vector(sub_sample %*% gamma[1:7])
#   a_grid <- seq(min(a), max(a), length.out = G_grid)
#   Qg <- plogis(-outer(a_grid, gamma[8] * y_sample, "+"))
# 
#   h <- approx(a_grid, as.vector(Qg %*% p), xout = a)$y
#   if (any(!is.finite(h)) || any(h <= 0)) return(1e10)
# 
#   term1 <- sum(log(Pi))
#   term2 <- sample_size * log(sum(Pi * p))
# 
#   -(term1 + wt*sum(log(h)))
# 
# }


gradient_neg_logLik_gamma <- function(gamma, Z_matrix, sub_sample, p, y_sample, sample_size, wt) {
  pq <- Pi_Q(gamma, Z_matrix, sub_sample, y_sample)

  Pi <- pq$Pi
  Q <- pq$Q
  W <- Q * (1 - Q)

  h <- as.vector(Q %*% p)
  h1 <- -as.vector(W %*% p)
  hy <- -as.vector(W %*% (p*y_sample))

  S <- sum(Pi*p)
  # pw <- p*Pi *(1-Pi)

  term1 <- as.vector(crossprod(Z_matrix, 1 - Pi))
  # term2 <- sample_size*as.vector(crossprod(Z_matrix, pw))/S
  term3 <- wt * c(as.vector(crossprod(sub_sample, h1/h)), sum(hy/h))

  -(term1 + term3)
}

# gradient_neg_logLik_gamma <- function(gamma, Z_matrix, sub_sample, p, y_sample, sample_size, wt) {
#   Pi <- as.vector(plogis(Z_matrix %*% gamma))
# 
#   a <- as.vector(sub_sample %*% gamma[1:7])
#   a_grid <- seq(min(a), max(a), length.out = G_grid)
#   Qg <- plogis(-outer(a_grid, gamma[8] * y_sample, "+"))
#   Wg <- Qg * (1 - Qg)
# 
#   h <- approx(a_grid, as.vector(Qg %*% p), xout = a)$y
#   h1 <- approx(a_grid, -as.vector(Wg %*% p), xout = a)$y
#   hy <- approx(a_grid, -as.vector(Wg %*% (p*y_sample)), xout = a)$y
# 
# 
#   term1 <- as.vector(crossprod(Z_matrix, 1 - Pi))
#   term3 <- wt * c(as.vector(crossprod(sub_sample, h1/h)), sum(hy/h))
# 
#   -(term1 + term3)
# }
# 
# prep_p <- function(gamma, sub_sample, y_sample) {
#   a <- as.vector(sub_sample %*% gamma[1:7])
#   N <- length(a)
#   a_grid <- seq(min(a), max(a), length.out = G_grid)
#   pos <- (a - a_grid[1])/(a_grid[2] - a_grid[1]) + 1
#   lo <- as.integer(pmin(floor(pos), G_grid - 1L))
#   t <- pos - lo
# 
#   A <- sparseMatrix(i = c(lo, lo + 1L), j = c(seq_len(N), seq_len(N)),
#                     x = c(1 - t, t), dims = c(G_grid, N))
#   Qg <- plogis(-outer(a_grid, gamma[8] * y_sample, "+"))
#   list(A = A, Qg = Qg)
# }
# 
# h_units <- function(p, prep) as.vector(crossprod(prep$A, prep$Qg %*% p))

## Likelihood for p_i optimization, where Gamma is assumed to be known
neg_logLik_p <- function(p, Q, wt) {

  if (any(p <= 0)) return(1e10)

  h <- as.vector(Q %*% p)

  term1 <- sum(log(p))
  # term2 <- sample_size * log(sum(p*Pi))
  -(term1 + wt*sum(log(h)))
}

# neg_logLik_p <- function(p, prep, sample_size, wt) {
#   if (any(p <= 0)) return(1e10)
#   -(sum(log(p)) + sum(log(h_units(p, prep))))
# }

gradient_neg_logLik_p <- function(p, Q, wt) {

  h <- as.vector(Q %*% p)
  # S <- sum(p*Pi)
  # term2 <- sample_size*Pi/S

  -(1/p + wt*as.vector(crossprod(Q, 1/h)))
}

# gradient_neg_logLik_p <- function(p, prep, sample_size, wt) {
#   h <- h_units(p, prep)
#   v <- as.vector(prep$A %*% (1 / h))
#   -(1 / p + as.vector(crossprod(prep$Qg, v)))
# }

## Equality constraints
heq <- function(p, ...) {
  c(
    sum(p) - 1,
    sum(p * final_sample_pref_raw$NDVI) - mean_pop_NDVI,
    sum(p * final_sample_pref_raw$EVI) - mean_pop_EVI,
    sum(p * final_sample_pref_raw$DEM) - mean_pop_DEM,
    sum(p * final_sample_pref_raw$SLOPE) - mean_pop_SLOPE,
    sum(p * final_sample_pref_raw$x) -  mean_pop_x,
    sum(p * final_sample_pref_raw$y) - mean_pop_y
  )
}

heq.jac <- function(p, ...) {
  rbind(
    rep(1, length(p)),
    final_sample_pref_raw$NDVI,
    final_sample_pref_raw$EVI,
    final_sample_pref_raw$DEM,
    final_sample_pref_raw$SLOPE,
    final_sample_pref_raw$x,
    final_sample_pref_raw$y
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
set.seed(123)
gamma_init <- rnorm(ncol(Z_matrix), mean = 0, sd = 0.1)
p_init <- rep(1/sample_size, sample_size)


tol_mu <- 1e-3      ## Threshold for stopping when convergence is achieved
tol_LL <- 1e-3
tol_gamma <- 1e-4
tol_p <- 1e-5
max_iter <- 20

## Storing maximum change in p and gamma values, log-Likelihoods of gamma and p 
p_change <- numeric(max_iter)
gamma_change <- numeric(max_iter)
LL_gamma <- numeric(max_iter)
LL_p <- numeric(max_iter)
mu_PL <- numeric(max_iter)
LL_full <- numeric(max_iter)


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
c(status = st$status, t = st$t, uniform = 1 / sample_size)

p_init <- st$p

## Optimization loop
system.time(for (iter in 1:max_iter) {

  cat("iter", iter, "\n")
  print(gamma_init)
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
    method = "BFGS",
    control = list(maxit = 200, trace = 1)
  )

  gamma_new <- gamma_fit$par
  print(gamma_new)
  
  cat("Gamma convergence:", gamma_fit$convergence, "LL gamma:", gamma_fit$value, "\n")

  cat("iter", iter, "Gamma optim done", "\n")
  
  pq <- Pi_Q(gamma_new, Z_matrix, sub_sample, y_sample)
  # ## p Optimization
  # p_fit <- auglag(
  #   par = p_init,
  #   fn = neg_logLik_p,
  #   gr = gradient_neg_logLik_p,
  #   Q = pq$Q,
  #   Pi = pq$Pi,
  #   sample_size = sample_size,
  #   wt = wt,
  #   heq = heq,
  #   heq.jac = heq.jac,
  #   hin = hin,
  #   hin.jac = hin.jac
  # )

  # prep  <- prep_p(gamma_new, sub_sample, y_sample)
  # 
  # p_fit <- auglag(
  #   par = p_init, fn = neg_logLik_p, gr = gradient_neg_logLik_p,
  #   prep = prep, sample_size = sample_size, wt = wt,
  #   heq = heq, heq.jac = heq.jac, hin = hin, hin.jac = hin.jac
  # )
  
  p_fit <- newton_method(p_init, pq$Q, con$C, con$b, wt)
  p_new <- p_fit$p
  
  # p_new <- p_fit$par
  
  if (any(!is.finite(p_new))) {
    return(c(rep = rep_id, naive_mean = mean(final_sample_raw$AGB), mu_y_EL = NA))
  }
  
  p_new <- pmax(p_new, 1e-8)

  cat("iter", iter, "p_i optim done", "\n")

  

  cat("iter", iter, "| delta p =", delta_p, "| delta gamma =", delta_gamma, "\n")

  p_change[iter] <- delta_p
  gamma_change[iter] <- delta_gamma

  LL_p[iter] <- -p_fit$value
  LL_gamma[iter] <- -gamma_fit$value
  mu_PL[iter] <- sum(p_new * final_sample_pref_raw$AGB)
  LL_full[iter] <- full_logLik(gamma_new, p_new)

  
  if(iter > 1) {
    delta_mu <- abs(mu_PL[iter] - mu_PL[iter -1])
    delta_LL <- abs(LL_full[iter] - LL_full[iter - 1])
    delta_p <- max(abs(p_new - p_init))
    delta_gamma <- max(abs(gamma_new - gamma_init))
    if (delta_mu < tol_mu && delta_LL < tol_LL && delta_gamma < tol_gamma && delta_p < tol_p) {
      cat("Converged at iteration: ", iter, "\n")
      break
    }
  }

  gamma_init <- gamma_new
  p_init <- p_new
  

})

data.frame(iter = 1:max_iter, p_change, gamma_change, LL_p, LL_gamma, LL_full)

mu_PL

n_g <- length(gamma_new)
n_p <- length(p_new)

gamma_hat <- gamma_new
p_hat <- p_new
agb <- final_sample_pref_raw$AGB
mu_hat <- sum(agb* p_hat)



mu_hat_ll <- full_logLik(gamma_hat, p_hat)
cutoff <- mu_hat_ll - qchisq(0.95, 1)/2

make_constraints_mu <- function(p, mu, agb) {
  con <- make_constraints(p)
  m <- mean(agb)
  s <- sd(agb)
  list(C = rbind(con$C, (agb - m)/s), 
       b = c(con$b, (mu - m)/s))
}

start_p <- function(C, b, p) {
  r <- as.vector(C %*% p) - b
  p0 <- p - as.vector(crossprod(C, solve(tcrossprod(C), r)))
  if (all(p0 > 0)) return(p0)
  st <- interior_start(C, b)
  if (st$status == 0 && st$t > 1e-06) return(st$p)
}


fit_at_mu <- function(mu0, refit_gamma = TRUE, g = gamma_hat, p = p_hat, max_iter = 15, threshold = 1e-05) {
  con <- make_constraints_mu(p_hat, mu0, agb)
  p <- start_p(con$C, con$b, p_hat)
  if (is.null(p)) return(list(feasible = FALSE))
  
  g <- gamma_hat
  old_ll <- -Inf
  
  for (i in 1:max_iter) {
    
    if (refit_gamma) {
      ## Gamma optimization
      gamma_fit <- optim(
        par = g,
        fn = neg_logLik_gamma,
        gr = gradient_neg_logLik_gamma,
        Z_matrix = Z_matrix,
        sub_sample = sub_sample,
        y_sample = y_sample,
        p = p,
        sample_size = sample_size,
        wt = wt,
        method = "BFGS",
        control = list(maxit = 200)
      )
      g <- gamma_fit$par
    }
    

    pq <- Pi_Q(g, Z_matrix, sub_sample, y_sample)
    # ## p Optimization
    # p_fit <- auglag(
    #   par = p,
    #   fn = neg_logLik_p,
    #   gr = gradient_neg_logLik_p,
    #   Q = pq$Q,
    #   Pi = pq$Pi,
    #   sample_size = sample_size,
    #   wt = wt,
    #   heq = heq_mu,
    #   heq.jac = heq.jac_mu,
    #   hin = hin,
    #   hin.jac = hin.jac
    # )
    
    p_fit <- newton_method(p, pq$Q, con$C, con$b, wt)
    
    p <- p_fit$p
    
    ll <- full_logLik(g, p)
    if(!refit_gamma || abs(ll - old_ll) < threshold) break
    old_ll <- ll
    
  }
  
  list(iter = i, ll = ll, mean_bias = sum(p * agb) - mu0, feasible = TRUE)
}

f <- function(mu0, refit_gamma = TRUE) {
  r <- fit_at_mu(mu0, refit_gamma)
  
  if(!r$feasible) return(-100)
  r$ll - cutoff
}

system.time(r <- fit_at_mu(mu_hat))
c(f = r$ll - cutoff, bias = r$mean_bias, iter = r$iter)

test_mus <- mu_hat + c(-3, -2, -1, -0.5, 0.5, 1, 2, 3)
cbind(mu0 = test_mus, f = sapply(test_mus, f, refit_gamma = FALSE))

f_frozen <- function(mu) f(mu, refit_gamma = FALSE)
lower_fz <- uniroot(f_frozen, c(mu_hat - 2, mu_hat - 1),
                    f.lower = -2.2173613, f.upper = 0.8263768, tol = 0.001)$root
upper_fz <- uniroot(f_frozen, c(mu_hat + 1, mu_hat + 2),
                    f.lower = 0.6914238, f.upper = -3.2830270, tol = 0.001)$root

f_logging <- function(mu0) {
  t0 <- Sys.time()
  val <- f(mu0)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  line <- sprintf("[%s] mu0 = %.3f || f = %+.4f || (%.0f s)",
                  format(Sys.time(), "%H:%M:%S"), mu0, val, secs)
  cat(line, "\n")
  cat(line, "\n", file = "profile_log_CI.txt", append = TRUE)
  val
}

f_full <- function(mu0) f(mu0, refit_gamma = TRUE)
lower_full <- uniroot(f_logging, c(mu_hat - 3, lower_fz), tol = 0.01, maxiter = 20)$root
upper_full <- uniroot(f_logging, c(upper_fz, mu_hat + 3), tol = 0.01, maxiter = 20)$root
c(lower_full, upper_full)


full_logLik <- function(optim_params, Z_matrix, sub_sample, y_sample, wt){
  gamma_optim <- optim_params[1:n_g]
  p_optim <- optim_params[(n_g + 1):(n_g + n_p)]
  pq <- Pi_Q(gamma_optim, Z_matrix, sub_sample, y_sample)
  h <- as.vector(pq$Q %*% p_optim)
  
  -(sum(log(pq$Pi)) + sum(log(p_optim)) + sum(log(h)))
}

full_gradient <- function(optim_params, Z_matrix, sub_sample, y_sample, wt){
  gamma_optim <- optim_params[1:n_g]
  p_optim <- optim_params[(n_g + 1):(n_g + n_p)]
  pq <- Pi_Q(gamma_optim, Z_matrix, sub_sample, y_sample)
  
  gamma_grad <- gradient_neg_logLik_gamma(gamma_optim, Z_matrix, sub_sample, p_optim, y_sample, sample_size, wt)
  p_grad <- gradient_neg_logLik_p(p_optim, pq$Q, pq$Pi, sample_size, wt)
  
  c(gamma_grad, p_grad)
}

k <- ncol(sub_sample)
Q <- plogis(-outer(as.vector(sub_sample %*% gamma_hat[1:k]), gamma_hat[k + 1] * y_sample, "+"))

pq <- Pi_Q(gamma_hat, Z_matrix, sub_sample, y_sample)
Q <- pq$Q
h <- as.vector(Q %*% p_hat)

Hpp <- (diag(1/p_new^2) + crossprod(Q/h))

Hgg <- jacobian(function(g) gradient_neg_logLik_gamma(g, Z_matrix, sub_sample, p_hat, y_sample, sample_size, wt), 
                x = gamma_hat, method = "simple")

# Hgp <- jacobian(function(p) gradient_neg_logLik_gamma(gamma_hat, Z_matrix, sub_sample, p, y_sample, sample_size, wt), 
#                 x = p_hat, method = "simple")
# Hpg <- t(Hgp)


# Hpg <- jacobian(function(g) {
#   Qg <- Pi_Q(g, Z_matrix, sub_sample, y_sample)$Q
#   gradient_neg_logLik_p(p_hat, Qg, wt)
# }, x = gamma_hat, method = "simple")
# 
# Hgp <- t(Hpg)

Hpg <- jacobian(function(g) gradient_neg_logLik_p(p_new, Pi_Q(g, Z_matrix, sub_sample, y_sample)$Q, wt),
                x = gamma_hat, method = "simple")

Hgp <- t(Hpg)


H_block <- rbind(cbind(Hgg, Hgp), cbind(Hpg, Hpp))

H_block <- (H_block + t(H_block))/2

# H_new <- hessian(fun = full_logLik, x = optim_params, Z_matrix = Z_matrix, 
#                  sub_sample = sub_sample, y_sample = y_sample, wt = wt)


# H_new <- jacobian(optim_params, func = full_gradient, Z_matrix = Z_matrix, sub_sample = sub_sample, 
#                y_sample = y_sample, wt = wt, method = "simple")

# H <- optimHess(optim_params, fn = full_logLik, gr = full_gradient, Z_matrix = Z_matrix, sub_sample = sub_sample, 
#                y_sample = y_sample, wt = wt)

# C <- cbind(matrix(0, nrow(heq.jac(p_new)), n_g), heq.jac(p_new))

Cp <- heq.jac(p_hat)
Cs <- rbind(Cp[1, ], t(scale(t(Cp[-1, ]))))

C  <- cbind(matrix(0, nrow(Cs), n_g), Cs)

K <- rbind(cbind(H_block, t(C)),
           cbind(C, matrix(0, nrow(C), nrow(C))))

V <- solve(K)[1:(n_g + n_p), 1: (n_g + n_p)]

g_vec <- c(rep(0, n_g), final_sample_pref_raw$AGB)
mu_hat <- sum(p_hat*final_sample_pref_raw$AGB)
var_mu <- drop(t(g_vec) %*% V %*% g_vec)

mu_hat + 1.96 * sqrt(var_mu)
mu_hat - 1.96 * sqrt(var_mu)





final_pref_test <- final_sample_pref_raw
final_pref_test$p_fullLik <- p_fullLik
final_pref_test$p_oldLik <- p_new
final_pref_test$p_plDiff <- final_pref_test$p_oldLik - final_pref_test$p_fullLik
final_pref_test$p_HierLik <- p_hat
#################################################################################


