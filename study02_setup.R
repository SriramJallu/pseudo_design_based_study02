# Author: Sriram Jallu

library(terra)
library(sp)
library(spcosa)
library(future.apply)

paths <- list(
  data = "../study02_data",
  results = "../study02_results/study02_agg"
)