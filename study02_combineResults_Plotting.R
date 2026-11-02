library(ggplot2)

agb <- rast("../study02_results/study02_agg/study02_agb_simulated_agg_exp_NDVI_EVI_AGB_expCorrResid180.tif")
results_EL <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_acessSample_EL_var.csv")
results_NewLik <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_acessSample_NewLik_var.csv")


head(results_EL)
colnames(results_EL)[colnames(results_EL) == "ci_l"] <- "ci_l_EL"
colnames(results_EL)[colnames(results_EL) == "ci_u"] <- "ci_u_EL"

results_EL$ci_l_pseudoWT <- results_EL$pseudo_weight_mean - 1.96 * sqrt(results_EL$var_g)
results_EL$ci_u_pseudoWT <- results_EL$pseudo_weight_mean + 1.96 * sqrt(results_EL$var_g)


true_mean <- global(agb, "mean", na.rm = TRUE)[1, 1]

results_EL$coverage_pseudoWT <- results_EL$ci_l_pseudoWT <= true_mean & results_EL$ci_u_pseudoWT >= true_mean
results_EL$coverage_EL <- results_EL$ci_l_EL <= true_mean & results_EL$ci_u_EL >= true_mean
head(results_EL)


head(results_NewLik)
colnames(results_NewLik)[colnames(results_NewLik) == "mu_y_EL"] <- "mu_y_NewLik"
colnames(results_NewLik)[colnames(results_NewLik) == "var_mu"] <- "var_mu_NewLik"
colnames(results_NewLik)[colnames(results_NewLik) == "ci_l"] <- "ci_l_NewLik"
colnames(results_NewLik)[colnames(results_NewLik) == "ci_u"] <- "ci_u_NewLik"

colnames(results_NewLik)[colnames(results_NewLik) == "coverage"] <- "coverage_NewLik"
results_NewLik$coverage_NewLik <- results_NewLik$ci_l_NewLik <= true_mean & results_NewLik$ci_u_NewLik >= true_mean
head(results_NewLik)

results <- cbind(results_EL, results_NewLik[, !(names(results_NewLik) %in% c("rep", "naive_mean"))])

write.csv(results, "C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_acessSample_var.csv", row.names = FALSE)


clustSample_results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_clustSample_var.csv")
accessSample_results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_acessSample_var.csv")
prefSample05_results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_prefSample05_var.csv")
prefSample10_results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_prefSample10_var.csv")
prefSampleCovars10_results <- read.csv("C:/Documents_PhD/Study02/study02_results/study02_agg/study02_results_agbSimulation_exp_NDVI_EVI_AGB_expCorrResid180_prefSampleCovars10_var.csv")


## Accuracy metrics over the MC simulations
accuracy_metrics <- function(df, pop_mean) {
  
  bias_inc <- mean(df$pseudo_weight_mean, na.rm = TRUE) - pop_mean
  sd_inc <- var(df$pseudo_weight_mean, na.rm = TRUE)
  rmse_inc <- sqrt(mean((df$pseudo_weight_mean - pop_mean)^2, na.rm = TRUE))
  coverage_inc <- mean(df$coverage_pseudoWT, na.rm = TRUE)* 100
  mse_inc <- (mean(df$var_g))
  
  bias_hier <- mean(df$mu_y_EL, na.rm = TRUE) - pop_mean
  sd_hier <- var(df$mu_y_EL, na.rm = TRUE)
  rmse_hier <- sqrt(mean((df$mu_y_EL - pop_mean)^2, na.rm = TRUE))
  coverage_hier <- mean(df$coverage_EL, na.rm = TRUE)* 100
  mse_hier <- (mean(df$var_mu))
  
  bias_NewLik <- mean(df$mu_y_NewLik , na.rm = TRUE) - pop_mean
  sd_NewLik <- var(df$mu_y_NewLik , na.rm = TRUE)
  rmse_NewLik <- sqrt(mean((df$mu_y_NewLik  - pop_mean)^2, na.rm = TRUE))
  coverage_NewLik <- mean(df$coverage_NewLik, na.rm = TRUE)* 100
  mse_NewLik <- (mean(df$var_mu_NewLik))
  
  return(data.frame(
    "Inc Bias" = round(bias_inc, 3),
    "Inc MC var" = round(sd_inc, 2),
    "Inc mean var" = round(mse_inc, 2),
    "Inc RMSE" = round(rmse_inc, 2),
    "Inc Coverage (%)" = round(coverage_inc, 2),
    
    "Hier Bias" = round(bias_hier, 3),
    "Hier MC var" = round(sd_hier, 2),
    "Hier mean var" = round(mse_hier, 2),
    "Hier RMSE" = round(rmse_hier, 2),
    "Hier Coverage (%)" = round(coverage_hier, 2),

    
    "NewLik Bias" = round(bias_NewLik, 3),
    "NewLik MC var" = round(sd_NewLik, 2),
    "NewLik mean var" = round(mse_NewLik, 2),
    "NewLik RMSE" = round(rmse_NewLik, 2),
    "NewLik Coverage (%)" = round(coverage_NewLik, 2),
    check.names = FALSE
  ))
}

clustSample_AM <- accuracy_metrics(clustSample_results, true_mean)
accessSample_AM <- accuracy_metrics(accessSample_results, true_mean)
prefSample05_AM <- accuracy_metrics(prefSample05_results, true_mean)
prefSample10_AM <- accuracy_metrics(prefSample10_results, true_mean)
prefSampleCovars10_AM <- accuracy_metrics(prefSampleCovars10_results, true_mean)


results_table <- rbind(
  cbind("Sample Type" = "Strongly Clustered", clustSample_AM),
  cbind("Sample Type" = "Accessibility Roads", accessSample_AM),
  cbind("Sample Type" = "Preferential 0.5", prefSample05_AM),
  cbind("Sample Type" = "Preferential 1.0", prefSample10_AM),
  cbind("Sample Type" = "Preferential 1.0, Covars", prefSampleCovars10_AM))

results_table


plot_df <- bind_rows(
  clustSample_results %>%
    select(naive_mean, mu_y_EL, mu_y_NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Strongly Clustered"),
  
  accessSample_results %>%
    select(naive_mean, mu_y_EL, mu_y_NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Accessibility Roads"),
  
  prefSample05_results %>%
    select(naive_mean, mu_y_EL, mu_y_NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Prefrential 0.5"),
  
  prefSample10_results %>%
    select(naive_mean, mu_y_EL, mu_y_NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Prefrential 1.0"),
  
  prefSampleCovars10_results %>%
    select(naive_mean, mu_y_EL, mu_y_NewLik) %>%
    pivot_longer(
      cols = everything(),
      names_to = "estimator",
      values_to = "estimate"
    ) %>%
    mutate(sampling = "Prefrential, with Covariates 1.0")
)

plot_df$estimator <- factor(plot_df$estimator, levels = c("naive_mean", "mu_y_EL", "mu_y_NewLik"))
est_labels <- c("naive_mean" = "Naive", "mu_y_EL" = "SBEL", "mu_y_NewLik" = "PIEL")

## Box plot and violin plot over the MC simulations
ggplot(plot_df, aes(x = estimator, y = estimate, fill = estimator, color = estimator)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_jitter(width = 0, alpha = 0.4, size = 2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  geom_hline(yintercept = true_mean, linetype = "dashed", color = "red") +
  facet_wrap(~sampling) +
  scale_x_discrete(labels = est_labels) +
  scale_fill_discrete(labels = est_labels) +
  scale_color_discrete(labels = est_labels) +
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

var_long <- function(df, design) {
  one <- function(est, v, cover, method) {
    mc <- round(var(est, na.rm = TRUE), 3)
    mv <- round(mean(v, na.rm = TRUE), 3)
    data.frame("Sampling Type" = design, Method = method,
               "MC Var" = mc, "Mean est. Var" = mv,
               Ratio = round(mv / mc, 3),
               "Coverage (95%)" = mean(cover, na.rm = TRUE) * 100,
               check.names = FALSE)
  }
  rbind(
    one(df$mu_y_EL, df$var_mu, df$coverage_EL, "SBEL"),
    one(df$mu_y_NewLik, df$var_mu_NewLik, df$coverage_NewLik, "PIEL")
  )
}

sets <- list("Strongly Clustered" = clustSample_results,
             "Accessibility Roads" = accessSample_results,
             "Preferential 0.5" = prefSample05_results,
             "Preferential 1.0" = prefSample10_results,
             "Preferential 1.0, Covars" = prefSampleCovars10_results)

var_table <- do.call(rbind, lapply(names(sets), function(nm) var_long(sets[[nm]], nm)))
var_table

ggplot(var_table, aes(Method, Ratio, fill = Method)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ "Sampling Type") +
  labs(y = "Mean estimated variance / MC variance", x = NULL) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
