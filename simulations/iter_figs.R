# Load simulations
# source("simulations/iter.R")


# ============================================================
# Plotting bias
# ============================================================
violin_n150 <- list()
violin_n150_filename <- list("violin_150_pse_bias.pdf", "violin_150_jnd_bias.pdf")

# bias

violin_n150[["pse_bias"]] <- ggplot(data = results, mapping = aes(y = bias_pse, x = method)) +
  geom_violin(quantile.linetype = TRUE, quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "PSE bias", x = NULL) +
  #coord_cartesian(ylim = c(-10, 10))+
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

violin_n150[["jnd_bias"]] <- ggplot(data = results, mapping = aes(y = bias_jnd, x = method)) +
  geom_violin(quantile.linetype = TRUE, quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "JND bias", x = NULL) +
  coord_cartesian(ylim = c(-3, 10)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")

# map2(.x = violin_n150_filename, .y = violin_n150, .f = ggsave, path = "Figs")


# par

violin_n150[["pse_fit"]] <- ggplot(data = results, mapping = aes(y = fit_pse, x = method)) +
  geom_violin(quantile.linetype = TRUE, quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "PSE estimates", x = NULL) +
  #coord_cartesian(ylim = c(-10, 10))+
  geom_hline(yintercept = true_pse, color = "red", linetype = "dashed")

violin_n150[["jnd_fit"]] <- ggplot(data = results, mapping = aes(y = fit_jnd, x = method)) +
  geom_violin(quantile.linetype = TRUE, quantiles = c(0.25, 0.5, 0.75)) +
  labs(y = "JND estimates", x = NULL) +
  coord_cartesian(ylim = c(0, 18)) +
  geom_hline(yintercept = true_jnd, color = "red", linetype = "dashed")

# map2(.x = violin_n150_filename, .y = violin_n150, .f = ggsave, path = "Figs")

# export

library(patchwork)

combined_plot_bias <- violin_n150[["jnd_bias"]] / violin_n150[["pse_bias"]]
combined_plot_bias <- combined_plot_bias + plot_annotation(tag_levels = 'A')

ggsave(filename = "combined_violins_bias.pdf", path = "Figs",
       plot = combined_plot_bias, device = "pdf", width = 6, height = 9)

# estimates 
combined_plot_fit <- violin_n150[["jnd_fit"]] / violin_n150[["pse_fit"]]
combined_plot_fit <- combined_plot_fit + plot_annotation(tag_levels = 'A')

ggsave(filename = "combined_violins_estimates.pdf", path = "Figs",
       plot = combined_plot_fit, device = "pdf", width = 6, height = 9)
