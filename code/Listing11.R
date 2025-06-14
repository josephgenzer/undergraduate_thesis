library(longmemo)
library(parallel)

n <- 2^16
h <- 2^12
H1 <- 0.9
H2 <- 0.9
centering1 <- 2 * H1 * log(h)
centering2 <- 2 * H2 * log(h)
theta <- pi / 3
P <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
R <- 2000
num_bins <- 20
num_cores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")

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

simulate_eigenvalues <- function(r) {
  Y <- simulate_Y(n, H1, H2, P)
  M_h <- calculate_M(Y, h)
  lambda <- sort(eigen(M_h, symmetric = TRUE)$values, decreasing = TRUE)
  transformed <- sqrt(n / h) * c(log(lambda[2]) - centering1,
                                 log(lambda[1]) - centering2)
  return(transformed)
}

eigen_matrix <- matrix(unlist(
  mclapply(seq_len(R), simulate_eigenvalues, mc.cores = num_cores, mc.preschedule = TRUE)
), ncol = 2, byrow = TRUE)

eigen_df <- data.frame(
  lambda_p_minus1 = eigen_matrix[, 1],
  lambda_p = eigen_matrix[, 2]
)

mu1 <- mean(eigen_df$lambda_p_minus1)
sigma1 <- sd(eigen_df$lambda_p_minus1)
mu2 <- mean(eigen_df$lambda_p)
sigma2 <- sd(eigen_df$lambda_p)

ks.test(eigen_df$lambda_p_minus1, "pnorm", mean = mu1, sd = sigma1)
ks.test(eigen_df$lambda_p, "pnorm", mean = mu2, sd = sigma2)

breaks1 <- seq(min(eigen_df$lambda_p_minus1), max(eigen_df$lambda_p_minus1), length.out = num_bins + 1)
breaks2 <- seq(min(eigen_df$lambda_p), max(eigen_df$lambda_p), length.out = num_bins + 1)

dens1 <- hist(eigen_df$lambda_p_minus1, breaks = breaks1, plot = FALSE, probability = TRUE)
dens2 <- hist(eigen_df$lambda_p, breaks = breaks2, plot = FALSE, probability = TRUE)

par(mfrow = c(1, 1), mar = c(5, 5.5, 2, 1), mgp = c(3.8, 0.7, 0))
hist(eigen_df$lambda_p_minus1,
     breaks = breaks1,
     probability = TRUE,
     main = "",
     xlab = expression(tilde(lambda)[p-1]),
     ylab = "",
     col = "#A3C1DA",
     border = "black",
     ylim = c(0, 1.1 * max(dens1$density)),
     axes = FALSE,
     cex.axis = 1.2,
     cex.lab = 1.3)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
mtext("Density", side = 2, line = 2.2, cex = 1.3)
box(bty = "o")

hist(eigen_df$lambda_p,
     breaks = breaks2,
     probability = TRUE,
     main = "",
     xlab = expression(tilde(lambda)[p]),
     ylab = "",
     col = "#A3C1DA",
     border = "black",
     ylim = c(0, 1.1 * max(dens2$density)),
     axes = FALSE,
     cex.axis = 1.2,
     cex.lab = 1.3)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
mtext("Density", side = 2, line = 2.2, cex = 1.3)
box(bty = "o")

eigen_df$norm_lambda_p_minus1 <- as.vector(scale(eigen_df$lambda_p_minus1))
eigen_df$norm_lambda_p <- as.vector(scale(eigen_df$lambda_p))

par(mar = c(6, 8, 2, 1), mgp = c(4, 0.7, 0))
plot(eigen_df$norm_lambda_p_minus1, eigen_df$norm_lambda_p,
     main = "",
     xlab = expression(tilde(lambda)[p-1] / sigma[p-1]),
     ylab = "",
     col = "#4682B4",
     pch = 16,
     cex = 0.6,
     xlim = c(-4, 4),
     ylim = c(-4, 4),
     axes = FALSE)
axis(1, at = seq(-4, 4, by = 2), cex.axis = 1)
axis(2, at = seq(-4, 4, by = 2), cex.axis = 1)
mtext(expression(tilde(lambda)[p] / sigma[p]), side = 2, line = 2, cex = 1.3)
box(bty = "o")

eigen_df$sum_squared <- eigen_df$norm_lambda_p^2 + eigen_df$norm_lambda_p_minus1^2
dens_sum <- hist(eigen_df$sum_squared, breaks = 30, plot = FALSE, probability = TRUE)

par(mar = c(6.5, 4, 2, 2), mgp = c(4.5, 0.7, 0))
hist(eigen_df$sum_squared,
     breaks = 30,
     probability = TRUE,
     main = "",
     xlab = expression((tilde(lambda)[p-1] / sigma[p-1])^2 + (tilde(lambda)[p] / sigma[p])^2),
     ylab = "",
     col = "#A3C1DA",
     border = "black",
     ylim = c(0, 1.1 * max(dens_sum$density)),
     cex.axis = 1.2,
     cex.lab = 1.0)
mtext("Density", side = 2, line = 2.2, cex = 1.3)
box(bty = "o")

ks.test(eigen_df$sum_squared, "pchisq", df = 2)

simulate_GOE_diff <- function(i) {
  M <- matrix(rnorm(4), nrow = 2)
  GOE <- (M + t(M)) / sqrt(2)
  eigs <- sort(eigen(GOE, symmetric = TRUE)$values, decreasing = TRUE)
  return(eigs[1] - eigs[2])
}

simulate_M_diff <- function(i) {
  Y <- simulate_Y(n, H1, H2, P)
  M_h <- calculate_M(Y, h)
  eigs <- sort(eigen(M_h, symmetric = TRUE)$values, decreasing = TRUE)
  diff <- eigs[1] - eigs[2]
  return(diff)
}

GOE_diffs <- unlist(mclapply(seq_len(R), simulate_GOE_diff, mc.cores = num_cores, mc.preschedule = TRUE))
Mh_diffs <- unlist(mclapply(seq_len(R), simulate_M_diff, mc.cores = num_cores, mc.preschedule = TRUE))

par(mgp = c(2.2, 0.7, 0))
qqplot(GOE_diffs, Mh_diffs,
       xlab = expression("2x2 GOE:" ~ lambda[2] - lambda[1]),
       ylab = expression("M(h):" ~ lambda[p] - lambda[p-1]),
       pch = 16,
       col = "#4682B4",
       cex.lab = 1.3,
       cex.axis = 1.2,
       axes = FALSE)
axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")

ks.test(scale(Mh_diffs), scale(GOE_diffs))
