# Author: Sriram Jallu

source("study02_setup.R")
source("study02_dataPrep.R")
source("study02_samplersScript.R")

files <- list(
  agb = file.path(paths$results, "study02_agb_simulated_agg_exp_NDVI_EVI_AGB_expCorrResid180.tif"),
  covars = file.path(paths$data, "study02_env_covariates#4_filtered.tif"),
  roads = file.path(paths$data, "study02_roadsClipped.shp")
)

config_settings <- list(
  sampler = "pref", ## pref, pref_covars, cluster
  sample_size = 500,
  alpha = 0,
  beta = 1,
  beta_covars = 1,
  cluster_strength = 0.9,
  one_per_pixel = TRUE,
  mc_reps = 50,
  seed = 123
)

data <- data_prepare(files$agb, files$covars, files$roads, agg_factor = 10, n_strata = 100, seed = 123)
packed_data <- pack_data(data)

set.seed(config_settings$seed)

test_sample <- sample_draw(config_settings, data)
head(test_sample)
