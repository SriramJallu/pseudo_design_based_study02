library(sp)
library(sf)
library(ggplot2)
library(terra)
library(spcosa)

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
mean(agb_df_plot$AGB)


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

