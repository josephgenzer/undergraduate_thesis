n <- 5000
p <- 2
R <- 2000
num_bins <- 20

simulate_eigenvalues <- function(n, R) {
  lambda2 <- numeric(R)
  for (i in 1:R) {
    X <- matrix(rnorm(n * p, mean = 0, sd = sqrt(c(1, 2))), nrow = n, ncol = p, byrow = TRUE)
    Sigma_hat <- (t(X) %*% X) / n
    lambda2[i] <- max(eigen(Sigma_hat, symmetric = TRUE)$values)
  }
  sqrt(n) * (lambda2 - 2)
}

asymp_lambda2 <- simulate_eigenvalues(n, R)

range <- range(asymp_lambda2)
breaks_seq <- seq(range[1], range[2], length.out = num_bins + 1)
densities <- hist(asymp_lambda2, breaks = breaks_seq, plot = FALSE, probability = TRUE)

hist(
  asymp_lambda2,
  breaks = breaks_seq,
  probability = TRUE,
  ylim = c(0, 1.1 * max(densities$density)),
  xlab = expression(tilde(lambda)[2]),
  ylab = "Density",
  main = "",
  col = "#A3C1DA",
  border = "black",
  axes = FALSE,
  cex.axis = 1.2,
  cex.lab = 1.3
)

axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")
