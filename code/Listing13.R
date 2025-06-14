library(MASS)
library(ggplot2)
library(longmemo)

n <- 2^16
h <- 2^12
H1 <- 0.9
H2 <- 0.9
theta <- pi / 3
P <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
cutoff <- 0.3
R <- 2000
beta <- 0.4

simulate_fBm <- function(n, H) {
  fGn <- simFGN0(n, H)
  fBm <- cumsum(fGn)
  return(fBm)
}

simulate_Y <- function(n, H1, H2, P) {
  X1 <- simulate_fBm(n, H1)
  X2 <- simulate_fBm(n, H2)
  X <- rbind(X1, X2)
  Y <- t(P %*% X)
  return(Y)
}

calculate_M <- function(Y, h) {
  n_h <- nrow(Y) - h
  differences <- Y[(h + 1):nrow(Y), ] - Y[seq_len(n_h), ]
  M <- crossprod(differences) / n_h
  return(M)
}

simulate_M_diff <- function(i) {
  Y <- simulate_Y(n, H1, H2, P)
  M_h <- calculate_M(Y, h)
  eigs <- sort(eigen(M_h, symmetric = TRUE)$values, decreasing = TRUE)
  diff <- eigs[1] - eigs[2]
  return(diff)
}

compute_beta_diff <- function(beta) {
  d <- rnorm(2, mean = 0, sd = sqrt(2))
  sub_diag <- sqrt(rchisq(1, df = beta))
  T_beta <- matrix(c(d[1], sub_diag, sub_diag, d[2]), nrow = 2)
  eigs <- sort(eigen(T_beta, symmetric = TRUE)$values, decreasing = TRUE)
  diff <- eigs[1] - eigs[2]
  return(diff)
}

eigdiff_M <- numeric(R)
eigdiff_beta <- numeric(R)

for (r in seq_len(R)) {
  eigdiff_M[r] <- simulate_M_diff(r)
  eigdiff_beta[r] <- compute_beta_diff(beta)
}

subset <- floor(cutoff * R)
eigdiff_M_subset <- sort(eigdiff_M)[seq_len(subset)]
eigdiff_beta_subset <- sort(eigdiff_beta)[seq_len(subset)]

eigdiff_M_norm <- scale(eigdiff_M_subset)
eigdiff_beta_norm <- scale(eigdiff_beta_subset)

par(mar = c(5, 5.5, 2, 1), mgp = c(3.1, 0.7, 0))
qqplot(
  eigdiff_beta_subset, eigdiff_M_subset,
  ylab = expression(M(h) ~ ":" ~ lambda[p] - lambda[p-1]),
  xlab = expression(paste(beta, "-Hermite matrix (", beta, " = 0.4): ", lambda[2] - lambda[1])),
  pch = 16,
  col = "#4682B4",
  cex.lab = 1.5,
  cex.axis = 1.2,
  axes = FALSE
)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")

KS_result <- ks.test(eigdiff_M_norm[, 1], eigdiff_beta_norm[, 1])
