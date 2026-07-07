library(sp)
library(sf)
library(ggplot2)
library(terra)
library(spcosa)
library(FNN)

setwd("C:/Documents_PhD/Study02/study02_data")


## Load the Target Variable
agb <- rast("AGB_Simulated.tif")
agb[agb < 0 | agb > 300] <- NA
plot(agb)
hist(agb)

## Prepare Spatial DF of the target variable
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
coordinates(agb_df) <- ~x+y
gridded(agb_df) <- TRUE

agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
sd(agb_df_plot$AGB)


## Strongly Clustered Samples done with Compact GeoStrata DOI: 10.1016/j.ecoinf.2022.101665
set.seed(123)
n_strata <- 100
stratification <- stratify(agb_df, nStrata = n_strata, nTry = 1)
plot(stratification)

strata_map <- as(stratification, "SpatialPixelsDataFrame")
strata_map_df <- as.data.frame(strata_map)


## Extracting stratas for clustering
strata_IDs <- unique(strata_map_df$stratumId)
cluster_strata <- sample(strata_IDs, 10)
other_strata <- setdiff(strata_IDs, cluster_strata)


## Drawing 90% sample size from the small group of stratas, and 10% from the large group of stratas
sample_size <- 500
n_cluster <- round(0.9*sample_size)
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
mean(final_sample$AGB)

## Plotting
final_sf <- st_as_sf(final_sample, coords = c("x", "y"), crs= crs(agb))

ggplot() +
  geom_raster(data = agb_df_plot, aes(x = x, y = y, fill = AGB)) +
  scale_fill_viridis_c(na.value = "transparent") +
  geom_sf(data = final_sf, aes(color = ID), size = 1.1, inherit.aes = FALSE) +
  scale_color_manual(
    name = "sample type",
    values = c("cluster_stratum" = "red", "other_stratum" ="black")
  ) +
  coord_sf() +
  theme_minimal()


## Reading Covariates
covars <- rast("Covariates.tif")

covars_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)
coordinates(covars_df) <- ~x+y
covars_names <- names(covars)


extracted <- extract(covars, final_sample[, c("x", "y")])
final_sample <- cbind(final_sample, extracted[, covars_names])
head(final_sample)

aux_df <- as.data.frame(covars, xy = TRUE, na.rm = TRUE)
head(aux_df)


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


## Creating a membership indicator for units in A
aux_norm_df$mem_indicator <- paste(aux_norm_df$x, aux_norm_df$y) %in% paste(final_sample_norm_df$x, final_sample_norm_df$y)
match_row <- match(paste(aux_norm_df$x, aux_norm_df$y), paste(final_sample_norm_df$x, final_sample_norm_df$y))


## Matrices for KNN
query_matrix <- as.matrix(aux_norm_df[, covars_names])
data_matrix <- as.matrix(final_sample_norm_df[, covars_names])


## KNN
k <- 2
knn_results <- get.knnx(data_matrix, query_matrix, k + 1) ## computes the ids of the closest K neighbors, along with distances, for units in query matrix
head(knn_results)


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
summary(sampling_model)




