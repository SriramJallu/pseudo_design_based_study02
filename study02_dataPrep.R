# Author: Sriram Jallu

data_prepare <- function(agb, covars, roads, agg_factor = 10, n_strata = 100, seed = 123) {
  
  agb <- rast(agb)
  covars <- rast(covars)
  covars_names <- names(covars)
  
  agb_df_pop <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
  
  agb_spdf <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
  coordinates(agb_spdf) <- ~x+y
  gridded(agb_spdf) <- TRUE
  
  
  valid_cells <- ifel(is.na(sum(covars)), NA, 1)
  covars <- aggregate(covars * valid_cells, fact = agg_factor, fun = "mean", na.rm = TRUE)
  names(covars) <- covars_names
  
  aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)
  
  agb_m <- global(agb, "mean", na.rm = TRUE)[1, 1]
  agb_s <- global(agb, "sd", na.rm = TRUE) [1, 1]
  agb_std <- (agb - agb_m)/agb_s
  
  covars_m <- global(covars, "mean", na.rm = TRUE)$mean
  covars_sd <- global(covars, "sd", na.rm = TRUE)$sd
  covars_std <- (covars - covars_m)/covars_sd
  
  roads <- sf::st_read(roads)
  roads <- sf::st_transform(roads, crs(covars))
  
  dist_to_roads <- distance(covars[[1]], vect(roads))
  names(dist_to_roads) <- "DistRoads"
  dist_to_roads <- mask(dist_to_roads, agb)
  
  d_mean <- global(dist_to_roads, "mean", na.rm = TRUE)[1, 1]
  d_sd <- global(dist_to_roads, "sd", na.rm = TRUE)[1, 1]
  dist_std <- (dist_to_roads - d_mean)/d_sd
  names(dist_std) <- "DistRoads_std"
  
  
  set.seed(seed)
  stratification <- stratify(agb_spdf, nStrata = n_strata, nTry = 1)
  strata_map <- as(stratification, "SpatialPixelsDataFrame")
  strata_map_df <- as.data.frame(strata_map)
  strata_IDs <- unique(strata_map_df$stratumId)
  
  list(agb = agb,
       agb_std = agb_std,
       covars = covars,
       covars_std = covars_std,
       covars_names = covars_names,
       dist_std = dist_std,
       aux_df = aux_df,
       strata_map_df = strata_map_df,
       strata_IDs = strata_IDs,
       true_mean = agb_m)
}

pack_data <- function(data) {
  lapply(data, function(x) if (inherits(x, "SpatRaster")) 
    terra::wrap(x) else x)
}

unpack_data <- function(data) {
  lapply(data, function(x) if (inherits(x, "PackedSpatRaster")) 
    terra::unwrap(x) else x)
}