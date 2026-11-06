# Author: Sriram Jallu

extract_agb <- function(xy, agb) {
  terra::extract(agb, xy)[, "AGB"]
}


poisson_points <- function(lambda, sample_size, agb, one_per_pixel = TRUE) {
  
  samp_intensity <- terra::values(lambda)
  samp_intensity <- samp_intensity * sample_size / sum(samp_intensity, na.rm = TRUE)
  
  valid <- !is.na(samp_intensity)
  counts <- integer(length(samp_intensity))
  counts[valid] <- rpois(sum(valid), samp_intensity[valid])
  
  if (one_per_pixel) counts <- as.integer(counts > 0)
  
  cells <- rep(seq_along(counts), times = counts)
  xy <- terra::xyFromCell(lambda, cells)
  
  data.frame(x = xy[, 1], y = xy[, 2], AGB = extract_agb(xy, agb))
}

pref_sample <- function(data, sample_size, alpha = 0, beta = 1, one_per_pixel = TRUE) {
  
  lambda <- exp(alpha + beta*data$agb_std)
  poisson_points(lambda, sample_size, data$agb, one_per_pixel)
}

pref_covars_sample <- function(data, sample_size, alpha = 0, beta_agb = 1, beta_covars = 1, covars = c("DEM", "SLOPE"),
                               one_per_pixel = TRUE) {
  
  lambda <- exp(alpha + beta_agb * data$agb_std - beta_covars * (sum(data$covars_std[[covars]])))
  poisson_points(lambda, sample_size, data$agb, one_per_pixel)
}

accessibility_sample <- function(data, sample_size, alpha = 0, beta = 1, one_per_pixel = TRUE) {
  
  lambda <- exp(alpha - beta*data$dist_std)
  poisson_points(lambda, sample_size, data$agb, one_per_pixel)
}

clustered_sample <- function(data, sample_size, strength = 0.9, n_cluster_strata = 10) {
  
  strata_df <- data$strata_map_df
  strata_IDs <- data$strata_IDs
  
  cluster_strata <- sample(strata_IDs, n_cluster_strata)
  
  n_cluster <- round(sample_size * strength)
  n_other <- sample_size - n_cluster
  
  cluster_df <- strata_df[strata_df$stratumId %in% cluster_strata, ]
  other_df <- strata_df[!(strata_df$stratumId %in% cluster_strata), ]
  
  cluster_pts <- cluster_df[sample.int(nrow(cluster_df), n_cluster), c("x", "y")]
  other_pts <- other_df[sample.int(nrow(other_df), n_other), c("x", "y")]
  
  xy <- rbind(cluster_pts, other_pts)
  data.frame(x = xy$x, y = xy$y, AGB = extract_agb(as.matrix(xy), data$agb))
}

sample_draw <- function(config, data) {
  switch(config$sampler,
         pref = pref_sample(data, config$sample_size, config$alpha, config$beta, one_per_pixel = config$one_per_pixel),
         pref_covars = pref_covars_sample(data, config$sample_size, config$alpha, beta_agb = config$beta,
                                          beta_covars = config$beta_covars, one_per_pixel =  config$one_per_pixel),
         accessibility = accessibility_sample(data, config$sample_size, config$alpha, config$beta,
                                              one_per_pixel =  config$one_per_pixel),
         cluster = clustered_sample(data, config$sample_size, strength =  config$cluster_strength),
         stop("Unknown sampler: ", config$sampler))
}

