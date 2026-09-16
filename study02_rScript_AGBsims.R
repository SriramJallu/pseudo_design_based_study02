library(sp)
library(sf)
library(ggplot2)
library(terra)
library(gstat)
library(moments)


setwd("C:/Documents_PhD/Study02/study02_data")


## Reading and Normalizing Covariates
covars <- rast("study02_env_covariates#4.tif")
agb <- rast("study02_agb_simulated_exp_NDVI_EVI.tif")
agb[agb == -9999] <- NA
print(covars)

common_mask <- !is.na(agb) & !any(is.na(covars))
agb_matched    <- mask(agb, common_mask, maskvalue = FALSE)
covars_matched <- mask(covars, common_mask, maskvalue = FALSE)


writeRaster(agb_matched, "study02_agb_CCI_filtered.tif")
writeRaster(covars_matched, "study02_env_covariates#4_filtered.tif")

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
hist(agb, na.rm = TRUE)
global(agb, fun = "mean", na.rm = TRUE)


agb_df_plot <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
sd(agb_df_plot$AGB)
mean(agb_df_plot$AGB)
hist(agb_df_plot$AGB)
summary(agb_df_plot$AGB)
skewness(agb_df_plot$AGB)

# writeRaster(agb, "study02_agb_simulated_lm.tif", overwrite = TRUE)


## Creating a Gaussian Field to add to the AGB map
set.seed(123)
w <- rast(agb)
values(w) <- rnorm(ncell(w))

plot(w)

k <- focalMat(
  w, d = 1000, type = "Gauss"
)


w <- focal(
  w, w = k, fun = sum, na.rm = TRUE
)

summary(w)

w <- (w - global(w, "mean", na.rm = TRUE)[1, 1]) / global(w, "sd", na.rm = TRUE)[1, 1]
plot(w)

summary(w)

sigma_spatial <- 5      #Controls the strength of the variance of the field
w <- sigma_spatial * w
names(w) <- "SpatialResidual"
plot(w)

agb_new <- agb + w
names(agb_new) <- "AGB"
plot(agb_new)

agb_new_df <- as.data.frame(agb_new, xy = TRUE, na.rm = TRUE)
mean(agb_new_df$AGB)

writeRaster(agb_new, "study02_agb_simulated_lm_resid05.tif", overwrite = TRUE)


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



################################################################################

## Spatially correlated error simulation

# Prediction grid
pts <- as.data.frame(agb, xy = TRUE, na.rm = FALSE)
coordinates(pts) <- ~x + y


# Variogram model
vgm_model <- vgm(
  psill = 35,
  model = "Exp",
  range = 10000,
  nugget = 0
)

# Kriging and Simulating a spatially correlated error field
set.seed(123)
g <- gstat(formula = z ~1, locations = ~x + y, dummy = TRUE, beta = 0, model = vgm_model, nmax = 30)
sim <- predict(g, newdata = pts, nsim = 1)


# Filtering Error df to match ncells
agb_df <- as.data.frame(agb, xy = TRUE, na.rm = TRUE)
sim_df <- as.data.frame(sim)


filtered_sim_df <- merge(agb_df[, c("x", "y")], sim_df[, c("x","y","sim1")], by = c("x","y"), all.x = TRUE)
summary(filtered_sim_df)
filtered_sim_df_copy <- filtered_sim_df
filtered_sim_df_copy$sim1 <- filtered_sim_df_copy$sim1 - mean(filtered_sim_df_copy$sim1, na.rm = TRUE)
summary(filtered_sim_df_copy)
write.csv(filtered_sim_df_copy, "study02_results_expCorrelatedResid_psill35.csv", row.names = FALSE)


## Adding the error with agb df

## Reading the error csv and convertint to raster
sim_df <- read.csv("study02_results_expCorrelatedResid_psill80.csv")
sim_rast <- terra::rast(sim_df[, c("x", "y", "sim1")], type = "xyz") # Not used further

agb_new <- covars[[1]]
values(agb_new) <- NA
names(agb_new) <- "AGB"

agb_df_full <- as.data.frame(agb, xy = TRUE, na.rm = FALSE)
agb_df_full$AGB_Trend <- agb_df_full[[3]]
agb_df_full <- merge(agb_df_full, sim_df[, c("x", "y", "sim1")], by = c("x", "y"), all.x = TRUE)
agb_df_full$AGB <- agb_df_full$AGB_Trend + agb_df_full$sim1

cells <- cellFromXY(agb_new, as.matrix(agb_df_full[, c("x", "y")]))
values(agb_new)[cells] <- agb_df_full$AGB

writeRaster(agb_new, "study02_agb_simulated_exp_NDVI_EVI_expCorrResid80.tif", overwrite = TRUE)


## Checking for variation explaination by LM 
agb_trend <- rast("study02_agb_simulated_exp_NDVI_EVI.tif")
agb_trend_df <- as.data.frame(agb_trend, xy = TRUE, na.rm = FALSE)

var_agb_trend <- var(agb_trend_df[[3]], na.rm = TRUE)
var_resid <- var(sim_df[[3]], na.rm = TRUE)

r2 <- var_agb_trend/(var_agb_trend + var_resid)
r2
