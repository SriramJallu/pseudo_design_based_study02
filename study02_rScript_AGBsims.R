rm(list=ls())
library(sp)
library(sf)
library(ggplot2)
library(terra)
library(gstat)


setwd("C:/Documents_PhD/Study02/study02_data")


## Reading and Normalizing Covariates
covars <- rast("study02_env_covariates#4.tif")

normalize <- function(r) {
  r_min <- global(r, "min", na.rm = TRUE)[1, 1]
  r_max <- global(r, "max", na.rm = TRUE)[1, 1]
  
  return((r - r_min)/(r_max - r_min))
}


NDVI <- normalize(covars$NDVI)
EVI <- normalize(covars$EVI)
DEM <- normalize(covars$DEM)
SLOPE <- normalize(covars$SLOPE)


## Simulating AGB Map
# agb <- 10 + (20*(exp(NDVI))) + (15*(exp(EVI))) - (25*DEM) + (15*SLOPE)
agb <- 10 + 50*(NDVI) + 50*(EVI) - 20*(DEM) + 40*(SLOPE)
names(agb) <- "AGB"

agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
sd(agb_df_plot$AGB)
mean(agb_df_plot$AGB)
hist(agb_df_plot$AGB)
summary(agb_df_plot$AGB)

writeRaster(agb, "study02_agb_simulated_lm.tif", overwrite = TRUE)


## Creating a Gaussian Field to add to the AGB map
set.seed(123)
w <- rast(agb)
values(w) <- rnorm(ncell(w))

k <- focalMat(
  w, d = 1000, type = "Gauss"
)

w <- focal(
  w, w = k, fun = sum, na.rm = TRUE
)

w <- (w - global(w, "min", na.rm = TRUE)[1, 1]) / global(w, "sd", na.rm = TRUE)[1, 1]

sigma_spatial <- 15      #Controls the strength of the variance of the field
w <- sigma_spatial * w
names(w) <- "SpatialResidual"

agb_new <- agb + w
names(agb_new) <- "AGB"
plot(agb_new)

agb_new_df <- as.data.frame(agb_new, xy = TRUE, na.rm = TRUE)
mean(agb_new_df$AGB)

writeRaster(agb_new, "study02_agb_simulated_lm_resid15.tif", overwrite = TRUE)


agb <- rast("study02_agb_simulated_exp_NDVI_EVI.tif")
agb_05 <- rast("study02_agb_simulated_exp_NDVI_EVI_resid05.tif")
agb_10 <- rast("study02_agb_simulated_exp_NDVI_EVI_resid10.tif")
agb_15 <- rast("study02_agb_simulated_exp_NDVI_EVI_resid15.tif")

agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
agb_05_df <- as.data.frame(agb_05, xy = TRUE, na.rm = TRUE)
agb_10_df <- as.data.frame(agb_10, xy = TRUE, na.rm = TRUE)
agb_15_df <- as.data.frame(agb_15, xy = TRUE, na.rm = TRUE)

mean(agb_df$AGB)
mean(agb_05_df$AGB)
mean(agb_10_df$AGB)
mean(agb_15_df$AGB)


ggplot() +
  geom_raster(data = agb_15_df, aes(x = x, y = y, fill = AGB)) +
  scale_fill_viridis_c(na.value = "transparent") +
  coord_sf() +
  theme_minimal()

