library(parallel)

n <- 5000
p <- 10
num_cores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")

compute_spectral_radius <- function(mat) {
  max(abs(eigen(mat, only.values = TRUE)$values))
}

spectral_radii <- numeric(n)

spectral_radii[2:n] <- unlist(mclapply(2:n, function(n_val) {
  total <- 0
  for (i in seq_len(n_val)) {
    X <- matrix(rnorm(n_val * p), nrow = n_val)
    Sigma_hat <- cov(X)
    Delta_n <- Sigma_hat - diag(p)
    if (any(is.infinite(Delta_n)) || any(is.na(Delta_n))) next
    total <- total + compute_spectral_radius(Delta_n)
  }
  total / n_val
}, mc.cores = num_cores))

plot(2:n, spectral_radii[2:n],
     type = "l", col = "#4682B4", lwd = 1.5,
     xlab = "Sample Size (n)",
     ylab = expression("Spectral Radius of " * Delta[n]))
