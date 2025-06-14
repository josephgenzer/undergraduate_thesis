library(MASS)
library(ggplot2)
library(parallel)

n <- 5000
p <- 2
R <- 2000
RNGkind("L'Ecuyer-CMRG")

lambda_diff_SCM <- function() {
  X <- matrix(rnorm(n * p), nrow = n)
  Sigma_hat <- (1 / n) * t(X) %*% X
  lambda <- sort(eigen(Sigma_hat, symmetric = TRUE)$values)
  sqrt(n) * (lambda[2] - lambda[1])
}

lambda_diff_GOE <- function() {
  M <- matrix(rnorm(p^2), nrow = p)
  GOE <- sqrt(p) * 0.5 * (M + t(M))
  lambda <- sort(eigen(GOE, symmetric = TRUE)$values)
  lambda[2] - lambda[1]
}

num_cores <- detectCores() - 1

scm <- unlist(mclapply(seq_len(R), function(i) lambda_diff_SCM(), mc.cores = num_cores))
goe <- unlist(mclapply(seq_len(R), function(i) lambda_diff_GOE(), mc.cores = num_cores))

qqplot(goe, scm,
       xlab = expression(GOE[2](eta) * ":" ~ lambda[2] - lambda[1]),
       ylab = expression(hat(Sigma)[n] * ":" ~ lambda[2] - lambda[1]),
       pch = 16,
       col = "#4682B4",
       cex.lab = 1.5,
       cex.axis = 1.2,
       axes = FALSE)

axis(1, cex.axis = 1.2)
axis(2, cex.axis = 1.2)
box(bty = "o")
