library(MASS)
library(ggplot2)
library(parallel)

n <- 1000
p <- n / 2
c <- p / n
R <- 100
num_bins <- 17

top_eigenvalue <- function(X) {
  Sigma_hat <- (1 / n) * t(X) %*% X
  max(abs(eigen(Sigma_hat, only.values = TRUE)$values))
}


values <- unlist(mclapply(seq_len(R), function(i) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  top_eigenvalue(X)
}, mc.cores = num_cores))

scaled <- (values - mean(values)) / sd(values)

range <- range(scaled)
breaks_seq <- seq(range[1], range[2], length.out = num_bins + 1)
densities <- hist(scaled, breaks = breaks_seq, plot = FALSE, probability = TRUE)

par(mfrow = c(1, 1), mar = c(5, 5.5, 2, 1), mgp = c(3.1, 0.7, 0))
hist(
  scaled,
  breaks = breaks_seq,
  probability = TRUE,
  ylim = c(0, 1.1 * max(densities$density)),
  xlim = c(min(breaks_seq), max(breaks_seq)),
  xlab = expression("Normalized " * lambda[p]),
  ylab = "Density",
  main = "",
  col = "#A3C1DA",
  border = "black",
  axes = FALSE,
  cex.axis = 1.2,
  cex.lab = 1.5
)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")

par(mar = c(5, 5.5, 2, 1), mgp = c(3.1, 0.7, 0))
qqnorm(
  scaled,
  main = "",
  xlab = "N(0, 1) Sample",
  ylab = expression("Normalized" ~ lambda[p]),
  col = "#4682B4",
  pch = 16,
  axes = FALSE,
  cex.lab = 1.5,
  cex.axis = 1.2
)
qqline(scaled, col = "#4682B4", lwd = 2, lty = 2)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")
