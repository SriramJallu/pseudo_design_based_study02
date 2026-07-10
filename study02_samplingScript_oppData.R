rm(list=ls())
library(sp)
library(sf)
library(ggplot2)
library(terra)
library(spcosa)
library(FNN)
library(future.apply)
library(spatstat.geom)
library(spatstat.random)
library(alabama)


setwd("C:/Documents_PhD/Study02/study02_data")


## Load the Target Variable
agb <- rast("study02_agb_simulated_exp_NDVI_EVI_corr.tif")
plot(agb)

## Prepare Spatial DF of the target variable
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
coordinates(agb_df) <- ~x+y
gridded(agb_df) <- TRUE

agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
sd(agb_df_plot$AGB)
mean(agb_df_plot$AGB)


## Standardize AGB values
m <- global(agb, "mean", na.rm = TRUE)[1, 1]
s <- global(agb, "sd", na.rm = TRUE) [1, 1]

agb_std <- (agb - m)/s


## Reading Covariates
covars <- rast("study02_env_covariates#4.tif")
covars_names <- names(covars)

aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)

## Wrap the Spatraster - for parallel computing later
covars_wrapped <- wrap(covars)
agb_wrapped <- wrap(agb)
agb_std_wrapped <- wrap(agb_std)


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
  final_sample$AGB <- extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}


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
  
  final_sample$AGB <- extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}

## Function for calculating pseudo weights using KNN imputation and modelling sampling mechanism (uses Logistic regression)
pseudo_weight_clac <- function(rep_id, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean) {
  
  library(FNN)
  library(terra)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)

  extracted <- extract(covars, final_sample[, c("x", "y")])
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
  
  
  ## Logistic regression of the membership indicator
  sampling_model <- glm(mem_indicator ~ NDVI + EVI + DEM + SLOPE + agb_impute, family = binomial(link = "probit"), data = aux_norm_df)
  
  non_prob_sample_pred <- final_sample_norm_df
  non_prob_sample_pred$agb_impute <- final_sample_norm_df$AGB
  
  
  ## Predicting inclusion probabilities & pseudo weights for the units of the non-probability sample
  final_sample_norm_df$inc_probs <- predict(sampling_model, newdata = non_prob_sample_pred, type = "response")

  final_sample_norm_df$pseudo_weights <- 1/final_sample_norm_df$inc_probs
  
  
  ## HT Estimator
  n <- nrow(aux_norm_df)

  naive_mean <- mean(final_sample_norm_df$AGB)
  pseudo_weight_mean <- (1/n) * sum(final_sample_norm_df$pseudo_weights * final_sample_norm_df$AGB)
  
  return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, final_sample_norm_df = final_sample_norm_df))
  
}


## Parallel computing: For repeating the estimation process over MC simulations
plan(multisession, workers = parallel::detectCores() - 1)

true_mean <- mean(agb_df_plot$AGB)
mc_reps <- 5

set.seed(123)
results_list <- future_lapply(
  1:mc_reps,
  function(r) {
    final_sample <- strongly_clustered_sample(s = 0.9, 500, agb_wrapped, strata_IDs)
    result <- pseudo_weight_clac(r, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
    cat(sprintf("MC iter %d done\n", r), file = "progress_log.txt", append = TRUE)
    result
    },
  future.seed = TRUE
)

plan(sequential)

results <- as.data.frame(do.call(rbind, results_list))
head(results)

## Accuracy metrics over the MC simulations
bias_weighted <- mean(results$pseudo_weight_mean, na.rm = TRUE) - true_mean
sd_weighted <- sd(results$pseudo_weight_mean, na.rm = TRUE)
rmse_weighted <- sqrt(mean((results$pseudo_weight_mean - true_mean)^2, na.rm = TRUE))

print(bias_weighted)
print(sd_weighted)
print(rmse_weighted)


################################################################################################################

## Set up for Empirical Likelihood, to estimate point masses (pseudo weights)
set.seed(123)

## Strongly Clustered Samples done with Compact GeoStrata DOI: 10.1016/j.ecoinf.2022.101665
set.seed(123)
n_strata <- 100
stratification <- stratify(agb_df, nStrata = n_strata, nTry = 1)
plot(stratification)

strata_map <- as(stratification, "SpatialPixelsDataFrame")
strata_map_df <- as.data.frame(strata_map)

## Extracting stratas for clustering
strata_IDs <- unique(strata_map_df$stratumId)

## Drawing an opportunistic sample
# final_sample_pref <- pref_sample(alpha = 0, beta = 1.2, sample_size = 500, agb_wrapped, agb_std_wrapped)
final_sample_pref <- strongly_clustered_sample(s = 0.9, sample_size = 500, agb_wrapped, strata_IDs)
final_sample_pref_sf <- st_as_sf(final_sample_pref, coords = c("x", "y"), crs = crs(agb))

## Calling the function for estimating weights - modelling sampling mechanism and KNN
pref_weights <- pseudo_weight_clac(1, final_sample_pref, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
pref_sample_udpate <- pref_weights$final_sample_norm_df

summary(lm(AGB ~ NDVI + EVI + DEM + SLOPE, data = pref_sample_udpate))$r.squared

## Regressing pseudo weights for smooth values
expected_weights_model <- glm(pseudo_weights ~ NDVI + EVI + DEM + SLOPE + AGB, family = Gamma(link = "log"), data = pref_sample_udpate)
summary(expected_weights_model)

pref_sample_udpate$expected_pseudo_weights <- predict(expected_weights_model, type = "response")


## Empirical Likelihood
sample_size <- nrow(pref_sample_udpate)
pseudo_weights <- pref_sample_udpate$pseudo_weights
expected_pseudo_weights <- pref_sample_udpate$expected_pseudo_weights

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

## Equality constraints
heq <- function(p) {
  c(
    sum(p) - 1,
    sum(p * pref_sample_udpate$NDVI_raw) - mean_NDVI_pop,
    sum(p * pref_sample_udpate$EVI_raw) - mean_EVI_pop,
    sum(p * pref_sample_udpate$DEM_raw) - mean_DEM_pop,
    sum(p * pref_sample_udpate$SLOPE_raw) - mean_SLOPE_pop
  )
}


## Inequality constraints
hin <- function(p) {
  p - 1e-8
}


## Initializing point masses
p_init <- rep(1/sample_size, sample_size)


## Augmented Lagrangian - for constrained optimization
result_EL <- auglag(par = p_init, fn = neg_logLik, heq = heq, hin = hin, 
                    control.outer = list(trace = FALSE))

p_hat <- result_EL$par  #point masses


## HT estimator using the weights
mu_y_EL <- sum(p_hat * pref_sample_udpate$AGB)
mu_y_EL

(mean(pref_sample_udpate$AGB))
(mean(agb_df_plot$AGB))

################################################################################################################

## Ploting the samples
ggplot() +
  geom_raster(data = agb_df_plot, aes(x = x, y = y, fill = AGB)) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = final_sample_cluster_sf, aes(color = ID), size = 1.1, inherit.aes = FALSE) +
  scale_color_manual(
    name = "sample type",
    values = c("cluster_stratum" = "red", "other_stratum" ="black")
  ) +
  coord_sf() +
  theme_minimal()


ggplot() +
  geom_raster(data = agb_df_plot, aes(x = x, y = y, fill = AGB)) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = final_sample_pref_sf, color = "red", size = 1.1, inherit.aes = FALSE) +
  coord_sf() +
  theme_minimal()






