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
agb <- 10 + (35*(exp(NDVI))) + (35*(exp(EVI))) - (20*DEM) + (40*SLOPE)


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

sigma_spatial <- 5      #Controls the strength of the variance of the field
w <- sigma_spatial * w
names(w) <- "SpatialResidual"

agb_new <- agb + w
names(agb_new) <- "AGB"

writeRaster(agb_new, "study02_agb_simulated_exp_NDVI_EVI_corr.tif", overwrite = TRUE)


