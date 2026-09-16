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
# agb <- rast("study02_agb_CCI_filtered.tif")
plot(agb)

agb_mean <- global(agb, fun = "mean", na.rm = TRUE)
print(agb_mean)

## Prepare Spatial DF of the target variable
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
coordinates(agb_df) <- ~x+y
gridded(agb_df) <- TRUE


agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
hist(agb_df_plot$AGB, main = "Population Distribution")


## Standardize AGB values
m <- global(agb, "mean", na.rm = TRUE)[1, 1]
s <- global(agb, "sd", na.rm = TRUE) [1, 1]

agb_std <- (agb - m)/s


## Reading Covariates
covars <- rast("../study02_data/study02_env_covariates#4.tif")
covars_names <- names(covars)
covars_names


aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)

## Wrap the Spatraster - for parallel processing
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
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  return(final_sample)
}

## Strongly Clustered Samples done with Compact GeoStrata DOI: 10.1016/j.ecoinf.2022.101665
set.seed(123)
n_strata <- 100
stratification <- stratify(agb_df, nStrata = n_strata, nTry = 1)
strata_map <- as(stratification, "SpatialPixelsDataFrame")
strata_map_df <- as.data.frame(strata_map)
strata_IDs <- unique(strata_map_df$stratumId)

strongly_clustered_sample_NM <- function(s = 0.9, sample_size, agb_wrapped, strata_IDs, strata_map_df) {
  
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

  naive_mean <- mean(final_sample$AGB)

  y_cluster <- mean(final_sample$AGB[final_sample$ID == "cluster_stratum"])
  y_other <- mean(final_sample$AGB[final_sample$ID == "other_stratum"])

  N_cluster <- nrow(cluster_df)
  N_other <- nrow(other_df)
  N_total <- N_cluster + N_other

  adjusted_mu <- (N_cluster/N_total)*y_cluster + (N_other/N_total)*y_other
  
  return(data.frame(naive_mean = naive_mean, adjusted_mu = adjusted_mu))
}


## Function for generating Preferential Samples - Strength is controlled by var - Beta
# pref_sample <- function(alpha = 0, beta = 1, sample_size, agb_wrapped, agb_std_wrapped) {
  
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


## Function for generating Conditioned Latin Hypercube Sampling
clhs_sample <- function(sample_size = 500, agb_wrapped, covars_df = covars_df_new){
  
  agb <- unwrap(agb_wrapped)
  
  clhs_result <- clhs(covars_df, size = sample_size, cost = "DistRoads", iter = 10000, simple = FALSE, progress = TRUE)
  
  sample_data <- clhs_result$sampled_data
  
  final_sample <- data.frame(sample_data[, c("x", "y")], row.names = NULL)
  
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


## Function for generating skewed samples - strength is controlled by alpha - not using at this moment
skewed_sample <- function(alpha = 5, sample_size, agb_wrapped, agb_std_wrapped) {
  
  agb <- unwrap(agb_wrapped)
  agb_std <- unwrap(agb_std_wrapped)
  
  v <- values(agb_std)[, 1]
  valid <- which(!is.na(v))
  v_valid <- v[valid]
  
  w <- pnorm(alpha * v_valid)
  w <- w/sum(w)
  
  cells <- sample(valid, sample_size, replace = FALSE, prob = w)
  
  xy <- xyFromCell(agb, cells)
  
  final_sample <- data.frame(x = xy[, 1], y = xy[, 2])
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  
  return(final_sample)
}


## Function for probability sample - Threshold based sampling on SLOPE and DEM
prob_sample <- function(sample_size = 100, agb_wrapped, agb_std_wrapped, covars_wrapped) {
  
  agb <- unwrap(agb_wrapped)
  agb_std <- unwrap(agb_std_wrapped)
  covars <- unwrap(covars_wrapped)
  
  dem_values <- values(covars$DEM)[, 1]
  slope_values <- values(covars$SLOPE)[, 1]
  agb_values <- values(agb_std)[, 1]
  
  dem_thresh <- quantile(dem_values, 0.3, na.rm = TRUE)
  slope_thresh <- quantile(slope_values, 0.3, na.rm = TRUE)
  
  access <- !is.na(agb_values) & !is.na(dem_values) & !is.na(slope_values) & (slope_values <= slope_thresh) & (dem_values <= dem_thresh)
  valid_cells <- which(access)
  
  cells <- sample(valid_cells, sample_size, replace = FALSE)
  
  xy <- xyFromCell(agb, cells)
  
  final_sample <- data.frame(x = xy[, 1], y = xy[, 2])
  final_sample$AGB <- terra::extract(agb, final_sample[, c("x", "y")])[, "AGB"]
  
  
  return(final_sample)
}


## Function for calculating pseudo weights using KNN imputation and modelling sampling mechanism (uses Logistic regression)
pseudo_weight_knn_clac <- function(rep_id, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean) {
  
  library(FNN)
  library(terra)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)

  extracted <- terra::extract(covars, final_sample[, c("x", "y")])
  final_sample <- cbind(final_sample, extracted[, covars_names])
  final_sample <- final_sample[complete.cases(final_sample[, covars_names]), ]##############
  
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
  # sampling_model <- glm(mem_indicator ~ NDVI + EVI + DEM + SLOPE + CHM + agb_impute, family = binomial(link = "logit"), data = aux_norm_df)
  # sampling_model <- gam(mem_indicator ~ s(NDVI) + s(EVI) + s(DEM) + s(SLOPE) + s(agb_impute), 
  #                       family = binomial(link = "logit"), data = aux_norm_df, method = "REML")
  
  non_prob_sample_pred <- final_sample_norm_df
  non_prob_sample_pred$agb_impute <- final_sample_norm_df$AGB
  
  
  ## Predicting inclusion probabilities & pseudo weights for the units of the non-probability sample
  final_sample_norm_df$inc_probs <- predict(sampling_model, newdata = non_prob_sample_pred, type = "response")
  # final_sample_norm_df$inc_probs <- predict(sampling_model, data = non_prob_sample_pred)$predictions[, 2]

  final_sample_norm_df$pseudo_weights <- 1/final_sample_norm_df$inc_probs
  
  
  ## HT Estimator
  n <- nrow(aux_norm_df)

  naive_mean <- mean(final_sample_norm_df$AGB)
  pseudo_weight_mean <- (1/n) * sum(final_sample_norm_df$pseudo_weights * final_sample_norm_df$AGB)
  
  # var_est <- VE.Hajek.Mean.NHT(VecY.s = final_sample_norm_df$AGB, VecPk.s = final_sample_norm_df$inc_probs, N = nrow(aux_norm_df))
  # std_est <- sqrt(var_est)
  # lower_95ci <- pseudo_weight_mean - (1.96 * std_est)
  # upper_95ci <- pseudo_weight_mean + (1.96 * std_est)
  # 
  # coverage <- (lower_95ci <= true_mean) & (upper_95ci >= true_mean)
  
  
  return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, final_sample_norm_df = final_sample_norm_df))
  # return(c(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, variance = var_est, std = std_est,
  #          lower_bound = lower_95ci, upper_bound = upper_95ci, coverage = coverage))
  
}

pseudo_weight_krige_clac <- function(rep_id, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean) {
  
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
  
  coordinates(final_sample_norm_df) <- ~x + y
  
  vg_agb <- gstat::variogram(AGB ~ 1, data = final_sample_norm_df)
  vgm_model_agb <- vgm(nugget = 150, psill = 150, range = 20000, model = "Sph")
  
  vgm_model_fit_agb <- fit.variogram(vg_agb, vgm_model_agb)
  
  
  aux_sp <- aux_norm_df
  coordinates(aux_sp) <- ~x + y
  kriged_agb <- krige(formula = AGB ~ 1, locations = final_sample_norm_df, newdata = aux_sp, model = vgm_model_fit_agb)
  
  aux_norm_df$agb_impute <- kriged_agb$var1.pred
  
  
  ## Logistic regression of the membership indicator
  sampling_model <- glm(mem_indicator ~ NDVI + EVI + DEM + SLOPE + agb_impute, family = binomial(link = "logit"), data = aux_norm_df)
  
  non_prob_sample_pred <- final_sample_norm_df
  non_prob_sample_pred$agb_impute <- final_sample_norm_df$AGB
  
  
  ## Predicting inclusion probabilities & pseudo weights for the units of the non-probability sample
  final_sample_norm_df$inc_probs <- predict(sampling_model, newdata = non_prob_sample_pred, type = "response")

  final_sample_norm_df$pseudo_weights <- 1/final_sample_norm_df$inc_probs
  
  
  ## HT Estimator
  n <- nrow(aux_norm_df)

  naive_mean <- mean(final_sample_norm_df$AGB)
  pseudo_weight_mean <- (1/n) * sum(final_sample_norm_df$pseudo_weights * final_sample_norm_df$AGB)
  
  var_est <- VE.Hajek.Mean.NHT(VecY.s = final_sample_norm_df$AGB, VecPk.s = final_sample_norm_df$inc_probs, N = nrow(aux_norm_df))
  std_est <- sqrt(var_est)
  lower_95ci <- pseudo_weight_mean - (1.96 * std_est)
  upper_95ci <- pseudo_weight_mean + (1.96 * std_est)
  
  coverage <- (lower_95ci <= true_mean) & (upper_95ci >= true_mean)
  
  
  # return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, final_sample_norm_df = final_sample_norm_df))
  return(c(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, variance = var_est, std = std_est,
           lower_bound = lower_95ci, upper_bound = upper_95ci, coverage = coverage))
  
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
  
  expected_weights_model <- glm(pseudo_weights ~ NDVI + EVI + DEM + SLOPE + AGB, family = Gamma(link = "log"), data = final_sample_norm_df)
  
  final_sample_norm_df$expected_pseudo_weights <- predict(expected_weights_model, type = "response")
  
  
  ## Empirical Likelihood
  sample_size <- nrow(final_sample_norm_df)
  pseudo_weights <- final_sample_norm_df$pseudo_weights
  expected_pseudo_weights <- final_sample_norm_df$expected_pseudo_weights
  
  mean_NDVI_pop <- mean(aux_df$NDVI)
  mean_EVI_pop <- mean(aux_df$EVI)
  mean_DEM_pop <- mean(aux_df$DEM)
  mean_SLOPE_pop <- mean(aux_df$SLOPE)
  
  ## Objective function of the EL
  neg_logLik <- function(p) {
    if (any(p <= 0)) return(1e10)
    term1 <- sample_size * log(sum(p*expected_pseudo_weights))
    term2 <- sum(log(p))
    -(term1 + term2)
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
  mu_y_EL <- sum(p_hat * final_sample_norm_df$AGB)
  
  
  # return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, mu_y_EL = mu_y_EL,
  #             final_sample_norm_df = final_sample_norm_df))
  return(c(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, mu_y_EL = mu_y_EL))
  
}


pseudo_weight_prob_knn_clac <- function(rep_id, final_sample, prob_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean) {
  
  library(FNN)
  library(terra)
  
  agb <- unwrap(agb_wrapped)
  covars <- unwrap(covars_wrapped)
  
  extracted <- terra::extract(covars, final_sample[, c("x", "y")])
  final_sample <- cbind(final_sample, extracted[, covars_names])
  final_sample <- final_sample[complete.cases(final_sample[, covars_names]), ]
  
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
  
  extracted_srs <- terra::extract(covars, prob_sample[, c("x", "y")])
  prob_sample <- cbind(prob_sample, extracted_srs[, covars_names])
  prob_sample_norm_df <- norm_func(prob_sample, covars_names, B_min, B_max)
  
  
  ## Creating a membership indicator for units in A
  aux_norm_df$mem_indicator <- paste(aux_norm_df$x, aux_norm_df$y) %in% paste(final_sample_norm_df$x, final_sample_norm_df$y)
  match_row <- match(paste(aux_norm_df$x, aux_norm_df$y), paste(final_sample_norm_df$x, final_sample_norm_df$y))
  
  
  donor_xy <- c(paste(final_sample_norm_df$x, final_sample_norm_df$y), paste(prob_sample_norm_df$x, prob_sample_norm_df$y))
  
  donor_agb <- c(final_sample_norm_df$AGB, prob_sample_norm_df$AGB)
  
  
  ## Matrices for KNN
  data_matrix <- as.matrix(rbind(final_sample_norm_df[, covars_names], prob_sample_norm_df[, covars_names]))
  
  match_row_donor <- match(paste(aux_norm_df$x, aux_norm_df$y), donor_xy)
  
  query_matrix <- as.matrix(aux_norm_df[, covars_names])
  
  
  ## KNN
  k <- 2
  knn_results <- get.knnx(data_matrix, query_matrix, k + 1) ## computes the ids of the closest K neighbors, along with distances, for units in query matrix
  
  
  ## Getting AGB values using the closest K neighbors, for B units, dropping its true y, so that there is a prediction there
  agb_impute <- numeric(nrow(aux_norm_df))
  
  for (i in seq_len(nrow(aux_norm_df))) {
    id <- knn_results$nn.index[i, ]
    
    if (!is.na(match_row_donor[i])) {
      id <- id[id != match_row_donor[i]][1:k]
    } else {
      id <- id[1:k]
    }
    
    agb_impute[i] <- mean(donor_agb[id])
  }
  
  
  aux_norm_df$agb_impute <- agb_impute
  
  s_row <- match(paste(prob_sample_norm_df$x, prob_sample_norm_df$y), paste(aux_norm_df$x, aux_norm_df$y))
  aux_norm_df$agb_impute[s_row] <- prob_sample_norm_df$AGB
  
  
  formula <- as.formula(paste("mem_indicator ~ x + y + agb_impute + ", paste(covars_names, collapse = " + ")))
  
  ## Logistic regression of the membership indicator
  sampling_model <- glm(formula, family = binomial(link = "logit"), data = aux_norm_df)
  # sampling_model <- glm(mem_indicator ~ NDVI + EVI + DEM + SLOPE + CHM + agb_impute, family = binomial(link = "logit"), data = aux_norm_df)
  # sampling_model <- gam(mem_indicator ~ s(NDVI) + s(EVI) + s(DEM) + s(SLOPE) + s(agb_impute), 
  #                       family = binomial(link = "logit"), data = aux_norm_df, method = "REML")
  
  non_prob_sample_pred <- final_sample_norm_df
  non_prob_sample_pred$agb_impute <- final_sample_norm_df$AGB
  
  
  ## Predicting inclusion probabilities & pseudo weights for the units of the non-probability sample
  final_sample_norm_df$inc_probs <- predict(sampling_model, newdata = non_prob_sample_pred, type = "response")
  # final_sample_norm_df$inc_probs <- predict(sampling_model, data = non_prob_sample_pred)$predictions[, 2]
  
  final_sample_norm_df$pseudo_weights <- 1/final_sample_norm_df$inc_probs
  
  
  ## HT Estimator
  n <- nrow(aux_norm_df)
  
  naive_mean <- mean(final_sample_norm_df$AGB)
  pseudo_weight_mean <- (1/n) * sum(final_sample_norm_df$pseudo_weights * final_sample_norm_df$AGB)
  
  # var_est <- VE.Hajek.Mean.NHT(VecY.s = final_sample_norm_df$AGB, VecPk.s = final_sample_norm_df$inc_probs, N = nrow(aux_norm_df))
  # std_est <- sqrt(var_est)
  # lower_95ci <- pseudo_weight_mean - (1.96 * std_est)
  # upper_95ci <- pseudo_weight_mean + (1.96 * std_est)
  # 
  # coverage <- (lower_95ci <= true_mean) & (upper_95ci >= true_mean)
  
  
  return(list(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, final_sample_norm_df = final_sample_norm_df))
  # return(c(rep = rep_id, naive_mean = naive_mean, pseudo_weight_mean = pseudo_weight_mean, variance = var_est, std = std_est,
  #          lower_bound = lower_95ci, upper_bound = upper_95ci, coverage = coverage))
  
}


## Parallel computing: For repeating the estimation process over MC simulations
plan(multisession, workers = parallel::detectCores() - 1)

true_mean <- mean(agb_df_plot$AGB)
mc_reps <- 50

set.seed(123)
results_list <- future_lapply(
  1:mc_reps,
  function(r) {
    # final_sample <- pref_sample(alpha = 0, beta = 0.5, sample_size = 500, agb_wrapped, agb_std_wrapped)
    final_sample <- strongly_clustered_sample(s = 0.9, 500, agb_wrapped, strata_IDs)
    # final_sample <- clhs_sample(sample_size = 500, agb_wrapped, covars_df = covars_df_new)
    # final_sample <- skewed_sample(alpha = 15, 500, agb_wrapped, agb_std_wrapped)
    result <- pseudo_weight_knn_EL_clac(r, final_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
    cat(sprintf("MC iter %d done\n", r), file = "progress_log.txt", append = TRUE)
    result
    },
  future.seed = TRUE
)

plan(sequential)

results <- as.data.frame(do.call(rbind, results_list))
write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid35_clhsSample_EL.csv", row.names = FALSE)
# write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid35_prefSample05_EL.csv", row.names = FALSE)
results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample10_EL_PL.csv")


## Accuracy metrics over the MC simulations
accuracy_metrics <- function(df, pop_mean) {
  bias_weighted <- mean(df$mu_y_EL, na.rm = TRUE) - pop_mean
  sd_weighted <- sd(df$mu_y_EL, na.rm = TRUE)
  rmse_weighted <- sqrt(mean((df$mu_y_EL - pop_mean)^2, na.rm = TRUE))
  
  return(data.frame(
    Bias = bias_weighted,
    Std = sd_weighted,
    RMSE = rmse_weighted
  ))
}

results_accuracy <- accuracy_metrics(results, true_mean)
print(results_accuracy)

## Box plots showing naive mean and pseudo-weighted mean
# plot_df <- data.frame(
#   estimator = rep(c("naive_mean", "pseudo_weight_mean"), each = nrow(results)),
#   estimate = c(results$naive_mean, results$pseudo_weight_mean)
# )
# 
# ggplot(plot_df, aes(x = estimator, y = estimate, fill = estimator)) +
#   geom_violin(alpha = 0.5, trim = FALSE) +
#   geom_boxplot(width = 0.15, outlier.shape = NA) +
#   geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
#   scale_x_discrete(labels = c("naive_mean" = "Naive", "pseudo_weight_mean" = "Pseudo-Weighted")) +
#   labs(title = "Naive vs Pseudo-Weighted Mean Estimator",
#        subtitle = paste("Dashed line = true mean:", round(true_mean, 2)),
#        y = "AGB estimate", x = NULL) +
#   theme_minimal() +
#   theme(legend.position = "none")


plot_df <- data.frame(
  estimator = rep(c("naive_mean", "pseudo_weight_mean", "mu_y_EL"), each = nrow(results)),
  estimate = c(results$naive_mean, results$pseudo_weight_mean, results$mu_y_EL)
)

ggplot(plot_df, aes(x = estimator, y = estimate, fill = estimator)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("naive_mean" = "Naive", "pseudo_weight_mean" = "Pseudo-Weighted", "mu_y_EL" = "EL Weighted")) +
  labs(title = "Naive, Pseudo-Weighted, EL Mean Estimator",
       subtitle = paste("Dashed line = true mean:", round(true_mean, 2)),
       y = "AGB estimate", x = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
)


################################################################################################################
## Plots to compare the results between different sampling configurations

exp_resid10_clusSamp <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_clusteredSample_EL.csv")
exp_resid10_prefSamp05 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample05_EL.csv")
exp_resid10_prefSamp10 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample10_EL.csv")
exp_resid10_prefSamp15 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample15_EL.csv")
exp_resid10_clhsSamp <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_clhsSample_EL.csv")

# lm_resid10_clusSamp <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_lm_expCorrResid10_clusteredSample.csv")
# lm_resid10_prefSamp15 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_lm_expCorrResid10_prefSample15.csv")


## Combined df for different sampling configurations (same underlying AGB)
plot_df <- bind_rows(
  exp_resid10_clusSamp %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Clustered"),

  exp_resid10_clhsSamp %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Accessibility cLHS"),
  
  exp_resid10_prefSamp05 %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Preferential 0.5"),
  
  exp_resid10_prefSamp10 %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Preferential 1.0")
)

plot_df$estimator <- factor(plot_df$estimator, levels = c("naive_mean", "pseudo_weight_mean", "mu_y_EL"))

## Box plot and violin plot over the MC simulations
ggplot(plot_df, aes(x = estimator, y = estimate, fill = sampling)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  facet_wrap(~sampling) +
  scale_x_discrete(labels = c("naive_mean" = "Naive", "mu_y_EL" = "EL Weighted", "pseudo_weight_mean" = "Inclusion Prob Weighted")) +
  labs(title = "Naive vs Pseudo-Weighted Mean Estimator",
       subtitle = paste("Dashed line = true mean:", round(true_mean, 2)), ## need to take of the "true_mean" variable
       y = "AGB estimate", x = NULL) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
)


## Need to update it manually based on which AGB simulation is being used
true_mean_lm = 54.23
true_mean_exp = 55.68


## Accuracy metrics
exp_resid10_clusSamp_AM <- accuracy_metrics(exp_resid10_clusSamp, pop_mean = true_mean)
exp_resid10_prefSamp05_AM <- accuracy_metrics(exp_resid10_prefSamp05, pop_mean = true_mean)
exp_resid10_prefSamp10_AM <- accuracy_metrics(exp_resid10_prefSamp10, pop_mean = true_mean)
exp_resid10_prefSamp15_AM <- accuracy_metrics(exp_resid10_prefSamp15, pop_mean = true_mean)


results_table <- rbind(
  cbind(Model = "Exp - Strongly Clustered", exp_resid10_clusSamp_AM),
  cbind(Model = "Exp - Preferential 0.5", exp_resid10_prefSamp05_AM),
  cbind(Model = "Exp - Preferential 1.0", exp_resid10_prefSamp10_AM),
  cbind(Model = "Exp - Preferential 1.5", exp_resid10_prefSamp15_AM))

print(results_table)
################################################################################################################

## Testing
## Set up for Empirical Likelihood, to estimate point masses (pseudo weights)
set.seed(123)


## Drawing an opportunistic sample
final_prob_sample <- prob_sample(500, agb_wrapped, agb_std_wrapped, covars_wrapped)
final_sample_pref <- pref_sample(alpha = 0, beta = 1, sample_size = 500, agb_wrapped, agb_std_wrapped)
# final_sample_pref <- strongly_clustered_sample(s = 0.90, sample_size = 500, agb_wrapped, strata_IDs)
# final_sample_pref <- clhs_sample(500, agb_wrapped, covars_df_new)
# final_sample_pref <- skewed_sample(alpha = 5, 500, agb_wrapped, agb_std_wrapped)
final_sample_pref_sf <- st_as_sf(final_sample_pref, coords = c("x", "y"), crs = crs(agb))

true_mean <- mean(agb_df_plot$AGB)
true_mean

## Calling the function for estimating weights - modelling sampling mechanism and KNN
pref_weights <- pseudo_weight_knn_clac(1, final_sample_pref, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
pref_weights <- pseudo_weight_prob_knn_clac(1, final_sample_pref, final_prob_sample, agb_wrapped, covars_wrapped, aux_df, covars_names, true_mean)
pref_sample_udpate <- pref_weights$final_sample_norm_df

# cor(pref_sample_udpate, method = "pearson")

summary(lm(AGB ~ NDVI + EVI + DEM + SLOPE, data = pref_sample_udpate))$r.squared

## Regressing pseudo weights for smooth values
expected_weights_model <- glm(pseudo_weights ~ NDVI + EVI + DEM + SLOPE + AGB, family = Gamma(link = "log"), data = pref_sample_udpate)
expected_weights_model <- gam(pseudo_weights ~ s(NDVI) + s(EVI) + s(DEM) + s(SLOPE) + s(CHM) + s(MSI) + s(AGB), 
                              family = Gamma(link = "log"), data = pref_sample_udpate, method = "REML")


pref_sample_udpate$expected_pseudo_weights <- predict(expected_weights_model, type = "response")

hist(pref_sample_udpate$expected_pseudo_weights, xlab = "1/Inclusion Probability", main = "")
hist(pref_sample_udpate$pseudo_weights)

## Empirical Likelihood
sample_size <- nrow(pref_sample_udpate)
pseudo_weights <- pref_sample_udpate$pseudo_weights
expected_pseudo_weights <- pref_sample_udpate$expected_pseudo_weights

mean_NDVI_pop <- mean(aux_df$NDVI)
mean_EVI_pop <- mean(aux_df$EVI)
mean_DEM_pop <- mean(aux_df$DEM)
mean_SLOPE_pop <- mean(aux_df$SLOPE)
# mean_MSI_pop <- mean(aux_df$MSI)
# mean_CHM_pop <- mean(aux_df$CHM)

## Objective function of the EL
neg_logLik <- function(p) {
  if (any(p <= 0)) return(1e10)
  term1 <- sample_size * log(sum(p*expected_pseudo_weights))
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

hist(p_hat)

## Mean estimator using the weights
mu_y_EL <- sum(p_hat * pref_sample_udpate$AGB)
mu_y_EL

(mean(pref_sample_udpate$AGB))
(mean(sum(pref_sample_udpate$AGB * pref_sample_udpate$pseudo_weights)))/nrow(aux_df)
(mean(agb_df_plot$AGB))
hist(pref_sample_udpate$AGB)
hist(agb_df_plot$AGB)

################################################################################################################

