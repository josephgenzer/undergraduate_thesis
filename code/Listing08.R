library(fracdiff)
library(ggplot2)
library(dplyr)
library(grid)
library(longmemo)

n <- 2^12
hurst_values <- c(0.2, 0.5, 0.8)

simulate_fbm <- function(n, H) {
  fgn <- simFGN0(n, H)
  cumsum(fgn)
}

data <- data.frame(
  time = rep(seq_len(n), times = length(hurst_values)),
  value = unlist(lapply(hurst_values, function(H) simulate_fbm(n, H))),
  hurst = factor(rep(hurst_values, each = n), levels = c(0.8, 0.5, 0.2))
)

levels(data$hurst) <- c("H = 0.8", "H = 0.5", "H = 0.2")

ggplot(data, aes(x = time, y = value)) +
  geom_line(size = 0.4, color = "#4682B4") +
  facet_wrap(~ hurst, ncol = 1, scales = "free") +
  labs(
    x = expression("Time " * (t)),
    y = expression(B[H](t))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, face = "plain", hjust = 0.5),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank()
  )
