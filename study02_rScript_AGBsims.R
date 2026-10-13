library(sp)
library(sf)
library(ggplot2)
library(terra)
library(gstat)
library(moments)
library(dplyr)


## Reading and Normalizing Covariates
covars <- rast("../study02_data/study02_env_covariates#4.tif")
agb <- rast("../study02_data/study02_agb_simulated_exp_NDVI_EVI.tif")
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
agb <- rast("../study02_data/study02_agb_simulated_exp_NDVI_EVI.tif")
valid_cells <- ifel(is.na(sum(covars)) | is.na(agb), NA, 1)

agb_agg <- aggregate(agb * valid_cells, fact = 10, fun = "mean", na.rm = TRUE)
names(agb_agg) <- "AGB"
agb <- agb_agg
writeRaster(agb, "../study02_data/study02_agg/study02_agb_simulated_exp_NDVI_EVI.tif")

agb <- rast("../study02_data/study02_agg/study02_agb_simulated_exp_NDVI_EVI.tif")
pts <- as.data.frame(agb, xy = TRUE, na.rm = FALSE)
coordinates(pts) <- ~x + y


# Variogram model
vgm_model <- vgm(
  psill = 30,
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
write.csv(filtered_sim_df_copy, 
  "../study02_results/study02_agg/study02_results_expCorrelatedResid_psill30.csv", row.names = FALSE)


## Adding the error with agb df

## Reading the error csv and convertint to raster
sim_df <- read.csv("../study02_results/study02_agg/study02_results_expCorrelatedResid_psill30.csv")
sim_rast <- terra::rast(sim_df[, c("x", "y", "sim1")], type = "xyz") # Not used further

agb_new <- agb
values(agb_new) <- NA
names(agb_new) <- "AGB"

agb_df_full <- as.data.frame(agb, xy = TRUE, na.rm = FALSE)
agb_df_full$AGB_Trend <- agb_df_full[[3]]
agb_df_full <- merge(agb_df_full, sim_df[, c("x", "y", "sim1")], by = c("x", "y"), all.x = TRUE)
agb_df_full$AGB <- agb_df_full$AGB_Trend + agb_df_full$sim1

cells <- cellFromXY(agb_new, as.matrix(agb_df_full[, c("x", "y")]))
values(agb_new)[cells] <- agb_df_full$AGB

writeRaster(agb_new, "../study02_results/study02_agg/study02_agb_simulated_agg_exp_NDVI_EVI_expCorrResid30.tif", overwrite = TRUE)


## Checking for variation explaination by LM 
agb_trend <- rast("../study02_data/study02_agg/study02_agb_simulated_exp_NDVI_EVI.tif")
agb_trend_df <- as.data.frame(agb_trend, xy = TRUE, na.rm = FALSE)

var_agb_trend <- var(agb_trend_df[[3]], na.rm = TRUE)
var_resid <- var(sim_df[[3]], na.rm = TRUE)

r2 <- var_agb_trend/(var_agb_trend + var_resid)
r2

agb_trend <- rast("../study02_results/study02_agg/study02_agb_simulated_agg_exp_NDVI_EVI_expCorrResid150.tif")
true_mean <- global(agb_trend, fun = "mean", na.rm = TRUE)
results1 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid30_prefSample10_EL.csv")
results2 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid30_prefSample10_EL_PL.csv")
results3 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid30_prefSample10_NewLik.csv")
# results3$mu_y_EL[which(results3$rep == 2)] <- 58.14673
# results3$mu_y_EL[which(results3$rep == 9)] <- 58.19389
# results3$mu_y_EL[which(results3$rep == 38)] <- 57.00733
# results3$mu_y_EL[which(results3$rep == 45)] <- 57.76768
# results3$mu_y_EL[which(results3$rep == 47)] <- 58.49563
# write.csv(results3, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample10_NewLik.csv", row.names = FALSE)

# results3$mu_y_EL[which(results3$rep == 7)] <- 56.28374
# write.csv(results3, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample05_NewLik.csv", row.names = FALSE)

# results3$mu_y_EL[which(results3$rep == 6)] <- 56.56093
# results3$mu_y_EL[which(results3$rep == 12)] <- 58.33056
# results3$mu_y_EL[which(results3$rep == 13)] <- 57.56930
# results3$mu_y_EL[which(results3$rep == 18)] <- 56.48056
# results3$mu_y_EL[which(results3$rep == 49)] <- 54.85079
# write.csv(results3, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSampleCovars10_NewLik.csv", row.names = FALSE)
# results2$mu_y_EL[which(results2$rep == 9)] <- 64.52978
# results2$mu_y_EL[which(results2$rep == 10)] <- 62.21529
# results2$mu_y_EL[which(results2$rep == 25)] <- 65.46388
# results2$mu_y_EL[which(results2$rep == 30)] <- 65.05274
# results2$mu_y_EL[which(results2$rep == 37)] <- 65.88242
# results2$mu_y_EL[which(results2$rep == 39)] <- 62.69499
# write.csv(results2, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSampleCovars10_EL_PL.csv", row.names = FALSE)


results2$mu_y_EL[which(results2$rep == 17)] <- 58.51281
write.csv(results2, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid30_prefSample10_EL_PL.csv", row.names = FALSE)



results <- results1
results$EL_PL <- results2$mu_y_EL
results$NewLik <- results3$mu_y_EL
write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid30_prefSample10.csv", row.names = FALSE)



results1 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample05.csv")
results2 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSample10.csv")
results3 <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_prefSampleCovars10.csv")
exp_resid10_clhsSamp <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_results_agbSimulation_exp_NDVI_EVI_expCorrResid150_clhsSample.csv")


## Combined df for different sampling configurations (same underlying AGB)
plot_df <- bind_rows(
  results1 %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean, EL_PL, NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Pref 0.5"),
  
  results2 %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean, EL_PL, NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Pref 1.0"),

  results3 %>%
    select(naive_mean, mu_y_EL, pseudo_weight_mean, EL_PL, NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Pref 1.0 Covars")
)

plot_df$estimator <- factor(plot_df$estimator, 
  levels = c("naive_mean", "pseudo_weight_mean", "mu_y_EL", "EL_PL", "NewLik"))

## Box plot and violin plot over the MC simulations
ggplot(plot_df, aes(x = estimator, y = estimate, fill = sampling)) +
  # geom_violin(alpha = 0.5, trim = FALSE) +
  geom_jitter(aes(color = estimator), width = 0, alpha = 0.4, size = 2) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  facet_wrap(~sampling) +
  scale_x_discrete(labels = c("naive_mean" = "Naive", "mu_y_EL" = "Hier Wt.", 
"pseudo_weight_mean" = "Inc Prob Wt.", "EL_PL" = "Prof Lik Wt.", "NewLik" = "Full Lik Wt.")) +
  labs(title = "Naive vs Pseudo-Weighted Mean Estimator",
       subtitle = paste("Dashed line = true mean:", round(true_mean, 2)), ## need to take of the "true_mean" variable
       y = "AGB estimate", x = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
)
