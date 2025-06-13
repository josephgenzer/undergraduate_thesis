library(MASS)
library(dplyr)
library(purrr)
library(ggplot2)
library(gridExtra)
library(grid)

n <- 100
c_vals <- c(1/10, 1/2, 1, 2)
num_bins <- 20

simulate_wishart <- function(n, p) {
  Sigma <- diag(p)
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S <- t(X) %*% X
  eigen(S, only.values = TRUE)$values
}

results <- list()

for (c in c_vals) {
  p <- floor(c * n)
  eigenvalues <- replicate(1000, simulate_wishart(n, p))
  eigenvalues_df <- data.frame(
    value = as.vector(eigenvalues),
    c = as.factor(rep(c, length(eigenvalues)))
  )
  results[[as.character(c)]] <- eigenvalues_df
}

combined_df <- do.call(rbind, results)

bin_info <- combined_df %>%
  group_by(c) %>%
  summarize(
    min_value = min(value),
    max_value = max(value),
    range_value = max_value - min_value,
    .groups = "drop"
  ) %>%
  mutate(
    margin = 0.1 * range_value,
    xlim_low = min_value - margin,
    xlim_high = max_value + margin,
    breaks_list = map2(min_value, max_value, ~ seq(.x, .y, length.out = num_bins + 1))
  )

combined_split <- split(combined_df, combined_df$c)

get_max_hist_density <- function(df, breaks) {
  h <- hist(df$value, breaks = breaks, plot = FALSE)
  max(h$density)
}

density_max_df <- map_df(names(combined_split), function(c_val) {
  df_sub <- combined_split[[c_val]]
  breaks_sub <- bin_info %>% filter(c == c_val) %>% pull(breaks_list) %>% .[[1]]
  ymax <- 1.1 * get_max_hist_density(df_sub, breaks_sub)
  data.frame(c = c_val, ymax = ymax)
})

plots <- lapply(names(combined_split), function(c_val) {
  df_sub <- combined_split[[c_val]]
  breaks_sub <- bin_info %>% filter(c == c_val) %>% pull(breaks_list) %>% .[[1]]
  ymax_sub <- density_max_df %>% filter(c == c_val) %>% pull(ymax)
  
  ggplot(df_sub, aes(x = value)) +
    geom_histogram(
      aes(y = after_stat(density)),
      breaks = breaks_sub,
      color = "black",
      fill = "#4682B4",
      alpha = 0.5,
      position = "identity"
    ) +
    coord_cartesian(ylim = c(0, ymax_sub)) +
    labs(x = NULL, y = NULL) +
    ggtitle(paste("c =", c_val)) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15, hjust = 0.5),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
    )
})

arranged_plots <- function(...) {
  plots <- list(...)
  grid.arrange(
    arrangeGrob(
      grobs = plots,
      ncol = 2,
      nrow = 2
    ),
    bottom = textGrob(
      expression("Eigenvalue (" * lambda * ")"),
      gp = gpar(fontsize = 1.5 * 11)
    ),
    left = textGrob(
      "Density", rot = 90,
      gp = gpar(fontsize = 1.5 * 11)
    )
  )
}

arranged_plots(plots[[1]], plots[[2]], plots[[3]], plots[[4]])
