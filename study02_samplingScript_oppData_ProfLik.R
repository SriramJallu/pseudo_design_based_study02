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


## Load the Target Variable
agb <- rast("../study02_data/study02_agb_simulated_exp_NDVI_EVI_expCorrResid150.tif")
plot(agb)

## Prepare Spatial DF of the target variable
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
coordinates(agb_df) <- ~x+y
gridded(agb_df) <- TRUE
agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)

## Standardize AGB values
m <- global(agb, "mean", na.rm = TRUE)[1, 1]
s <- global(agb, "sd", na.rm = TRUE) [1, 1]
agb_std <- (agb - m)/s

## Reading Covariates
covars <- rast("../study02_data/study02_env_covariates#4.tif")
covars_names <- names(covars)

aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)

## Wrap the Spatraster - for parallel processing
covars_wrapped <- wrap(covars)
agb_wrapped <- wrap(agb)
agb_std_wrapped <- wrap(agb_std)

## Function for generating Preferential Samples - Strength is controlled by var - Beta
pref_sample <- function(alpha = 0, beta = 1, sample_size, agb_wrapped, agb_std_wrapped) {
  
  agb <- unwrap(agb_wrapped)
  agb_std <- unwrap(agb_std_wrapped)
  
  ## Sampling Intensity - based on the Response
  lambda <- exp(alpha + beta * agb_std)
  v <- values(lambda)

  valid <- !is.na(v) & v > 0
  valid_cells <- which(valid)
  weights <- v[valid]
  
  cells <- sample(valid_cells, size = sample_size, replace = FALSE, prob = weights)
  
  ## Extracting XY and AGB values for the points
  xy <- xyFromCell(agb, cells)
  
  final_sample <- data.frame(
    x = xy[, 1],
    y = xy[, 2]
  )
  
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}


## Empirical Likelihood implementation, with Profile Likelihood
pesudo_weight_profileLik <- function(rep_id, final_sample, max_iter, agb_wrapped, covars_wrapped, aux_df, sample_size) {
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
  neg_logLik_gamma <- function(gamma, Z_matrix, p, sample_size){
    
    Pi <- P_i(gamma, Z_matrix)
    
    if (any(Pi <= 0 | Pi >= 1)) return(1e10)
    
    term1 <- sum(log(Pi))
    term2 <- sample_size * log(sum(Pi * p))
    
    -(term1 - term2)
  }
  
  ## Likelihood for p_i optimization, where Gamma is assumed to be known
  neg_logLik_p <- function(p, sample_size, gamma, Z_matrix) {
    
    if (any(p <= 0)) return(1e10)
    
    Pi <- P_i(gamma, Z_matrix)
    
    term1 <- sample_size * log(sum(p*Pi))
    term2 <- sum(log(p))
    -(term2 - term1)
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
  
  ## Inequality constraints
  hin <- function(p, ...) {
    p - 1e-8
  }
  
  ## Initialize gamma and p values for optimization
  gamma_init <- rnorm(ncol(Z_matrix), mean = 0, sd = 0.01)
  p_init <- rep(1/sample_size, sample_size)

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
      sample_size = sample_size,
      method = "BFGS"
    )
    
    ## Collecting the new gamma parameter values
    gamma_new <- gamma_fit$par

    ## p Optimization
    p_fit <- auglag(
      par = p_init,
      fn = neg_logLik_p,
      gamma = gamma_new,
      Z_matrix = Z_matrix,
      sample_size = sample_size,
      heq = heq,
      hin = hin,
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

## Parallel computing: For repeating the estimation process over MC simulations
plan(multisession, workers = parallel::detectCores() - 1)

true_mean <- mean(agb_df_plot$AGB)
mc_reps <- 25                       ## Monte Carlo simulations
sample_size <- 500

set.seed(123)
results_list <- future_lapply(
  1:mc_reps,
  function(r) {
    final_sample <- pref_sample(alpha = 0, beta = 1.5, sample_size = sample_size, agb_wrapped, agb_std_wrapped)
    result <- pesudo_weight_profileLik(r, final_sample, max_iter = 10, agb_wrapped, covars_wrapped, aux_df, sample_size)
    cat(sprintf("MC iter %d done\n", r), file = "progress_log.txt", append = TRUE)
    result
    },
  future.seed = TRUE
)

plan(sequential)

results <- as.data.frame(do.call(rbind, results_list))

write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample15_EL_PL.csv", row.names = FALSE)
results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample05_EL_PL.csv")


## Box plots
plot_df <- data.frame(
  estimator = rep(c("naive_mean", "mu_y_EL"), each = nrow(results)),
  estimate = c(results$naive_mean, results$mu_y_EL)
)

ggplot(plot_df, aes(x = estimator, y = estimate, fill = estimator)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("naive_mean" = "Naive", "mu_y_EL" = "EL Weighted")) +
  labs(title = "Naive, EL Mean Estimator",
       subtitle = paste("Dashed line = true mean:", round(true_mean, 2)),
       y = "AGB estimate", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none")


test <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample10_EL.csv")
summary(test)


#####################################################################
## Profile Likelihood Testing
set.seed(1234)
final_sample_pref <- pref_sample(0, 1.5, 500, agb_wrapped, agb_std_wrapped)
final_sample_pref_sf <- st_as_sf(final_sample_pref, coords = c("x", "y"), crs = crs(covars))
covar_vals <- extract(covars, final_sample_pref[, c("x", "y")])
final_sample_pref <- cbind(final_sample_pref, covar_vals[, -1])
final_sample_pref_raw <- final_sample_pref

final_sample_pref[colnames(final_sample_pref)] <- scale(final_sample_pref[colnames(final_sample_pref)])

Z_matrix <- as.matrix(cbind(1, final_sample_pref$x, final_sample_pref$y, final_sample_pref$NDVI, final_sample_pref$EVI, 
  final_sample_pref$DEM, final_sample_pref$SLOPE, final_sample_pref$AGB))

sample_size <- nrow(final_sample_pref)

## Initialize gamma and p values for optimization
gamma_init <- rnorm(ncol(Z_matrix), mean = 0, sd = 0.01)
p_init <- rep(1/sample_size, sample_size)

threshold <- 1e-06      ## Threshold for stopping when convergence is achieved
max_iter <- 15

## Storing maximum change in p and gamma values, log-Likelihoods of gamma and p 
p_change <- numeric(max_iter)
gamma_change <- numeric(max_iter)
LL_gamma <- numeric(max_iter)
LL_p <- numeric(max_iter)

## Optimization loop
set.seed(123)
for (iter in 1:max_iter) {

  cat("iter", iter, "\n")
  print(gamma_init)
  ## Gamma optimization
  gamma_fit <- optim(
    par = gamma_init,
    fn = neg_logLik_gamma,
    Z_matrix = Z_matrix,
    p = p_init,
    sample_size = sample_size,
    method = "BFGS"
  )

  gamma_new <- gamma_fit$par
  print(gamma_new)
  
  cat("Gamma convergence:", gamma_fit$convergence, "LL gamma:", gamma_fit$value, "\n")

  cat("iter", iter, "Gamma optim done", "\n")

  ## p Optimization
  p_fit <- auglag(
    par = p_init,
    fn = neg_logLik_p,
    gamma = gamma_new,
    Z_matrix = Z_matrix,
    sample_size = sample_size,
    heq = heq,
    hin = hin
  )

  p_new <- p_fit$par

  cat("iter", iter, "p_i optim done", "\n")

  delta_p <- max(abs(p_new - p_init))
  delta_gamma <- max(abs(gamma_new - gamma_init))

  cat("iter", iter, "| delta p =", delta_p, "| delta gamma =", delta_gamma, "\n")

  p_change[iter] <- delta_p
  gamma_change[iter] <- delta_gamma

  LL_p[iter] <- -p_fit$value
  LL_gamma[iter] <- -gamma_fit$value

  if (delta_p < threshold && delta_gamma < threshold) {
    cat("converged at iteration", iter, "\n")
    break
  }

  gamma_init <- gamma_new
  p_init <- p_new

}

range(p_new)
sum(p_new)
p_fit$convergence
p_fit$message

(heq(p_init))
summary(heq(p_new))

mu_y_EL_profile <- sum(p_new * final_sample_pref_raw$AGB)
mu_y_EL_profile
hist(LL_gamma)
hist(LL_p)
LL_p

LL_gamma35 <- LL_gamma
gamma_change35 <- gamma_change
LL_p35 <- LL_p
p_change35 <- p_change


