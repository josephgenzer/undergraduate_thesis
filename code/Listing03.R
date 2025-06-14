library(MASS)

n <- 250
p <- 10
Sigma <- diag(p)
num_bins <- 25

simulate_wishart <- function(n, p) {
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  S <- t(X) %*% X
  eigenvalues <- eigen(S, only.values = TRUE)$values
  return(eigenvalues)
}

eigenvalues <- as.vector(replicate(1000, simulate_wishart(n, p)))

range <- range(eigenvalues)
breaks_seq <- seq(range[1], range[2], length.out = num_bins + 1)

densities <- hist(
  eigenvalues,
  breaks = breaks_seq,
  plot = FALSE,
  probability = TRUE
)

hist(
  eigenvalues,
  breaks = breaks_seq,
  probability = TRUE,
  ylim = c(0, 1.1 * max(densities$density)),
  main = "",
  xlab = "Wishart Bulk Eigenvalues",
  ylab = "Density",
  col = "#A3C1DA",
  border = "black",
  axes = FALSE,
  cex.axis = 1.2,
  cex.lab = 1.3
)

axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")
