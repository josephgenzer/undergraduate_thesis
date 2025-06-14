library(fracdiff)
library(ggplot2)
library(gridExtra)
library(grid)

n <- 2^16
h_values <- c(2^4, 2^6, 2^8)
H1 <- 0.25
H2 <- 0.75
theta <- pi / 3
P <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)

simulate_fBm <- function(n, H) {
  fGn <- simFGN0(n, H)
  fBm <- cumsum(fGn)
  return(fBm)
}

simulate_Y <- function(n, H1, H2, P) {
  B1 <- simulate_fbm(n, H1)
  B2 <- simulate_fbm(n, H2)
  Y <- P %*% rbind(B1, B2)
  t(Y)
}

calculate_M <- function(Y, h) {
  n_h <- nrow(Y) - h
  diff <- Y[(h + 1):nrow(Y), ] - Y[seq_len(n_h), ]
  crossprod(diff) / n_h
}

estimate_alpha_diag <- function(H1, H2, h_values) {
  Y <- simulate_Y(n, H1, H2, P)
  log_h <- log(h_values)
  log_M11 <- numeric(length(h_values))
  log_M22 <- numeric(length(h_values))
  
  for (i in seq_along(h_values)) {
    M <- calculate_M(Y, h_values[i])
    log_M11[i] <- log(M[1, 1])
    log_M22[i] <- log(M[2, 2])
  }
  
  alpha_11 <- coef(lm(log_M11 ~ log_h))[2]
  alpha_22 <- coef(lm(log_M22 ~ log_h))[2]
  
  df <- data.frame(log_h, log_M11, log_M22)
  
  p1 <<- ggplot(df, aes(x = log_h, y = log_M11)) +
    geom_point(color = "black", size = 3) +
    geom_smooth(method = "lm", color = "#4682B4", se = FALSE, linewidth = 1.5) +
    labs(
      x = expression(log(h)),
      y = expression(log(M(h)[list(1,1)])),
      title = bquote(M(h)[list(1,1)] ~ ": " ~ alpha == .(sprintf("%.3f", alpha_11)))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    ylim(1, 9)
  
  p2 <<- ggplot(df, aes(x = log_h, y = log_M22)) +
    geom_point(color = "black", size = 3) +
    geom_smooth(method = "lm", color = "#4682B4", se = FALSE, linewidth = 1.5) +
    labs(
      x = expression(log(h)),
      y = expression(log(M(h)[list(2,2)])),
      title = bquote(M(h)[list(2,2)] ~ ": " ~ alpha == .(sprintf("%.3f", alpha_22)))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    ylim(1, 9)
  
  c(alpha_11 = alpha_11, alpha_22 = alpha_22)
}

estimate_alpha_eigen <- function(H1, H2, h_values) {
  Y <- simulate_Y(n, H1, H2, P)
  log_h <- log(h_values)
  log_lambda1 <- numeric(length(h_values))
  log_lambda2 <- numeric(length(h_values))
  
  for (i in seq_along(h_values)) {
    M <- calculate_M(Y, h_values[i])
    lambda <- sort(eigen(M)$values, decreasing = TRUE)
    log_lambda1[i] <- log(lambda[1])
    log_lambda2[i] <- log(lambda[2])
  }
  
  alpha_1 <- coef(lm(log_lambda1 ~ log_h))[2]
  alpha_2 <- coef(lm(log_lambda2 ~ log_h))[2]
  
  df <- data.frame(log_h, log_lambda1, log_lambda2)
  
  p3 <<- ggplot(df, aes(x = log_h, y = log_lambda1)) +
    geom_point(color = "black", size = 3) +
    geom_smooth(method = "lm", color = "#4682B4", se = FALSE, linewidth = 1.5) +
    labs(
      x = expression(log(h)),
      y = expression(log(lambda[2](M(h)))),
      title = bquote(lambda[2](M(h)) ~ ": " ~ alpha == .(sprintf("%.3f", alpha_1)))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    ylim(1, 9)
  
  p4 <<- ggplot(df, aes(x = log_h, y = log_lambda2)) +
    geom_point(color = "black", size = 3) +
    geom_smooth(method = "lm", color = "#4682B4", se = FALSE, linewidth = 1.5) +
    labs(
      x = expression(log(h)),
      y = expression(log(lambda[1](M(h)))),
      title = bquote(lambda[1](M(h)) ~ ": " ~ alpha == .(sprintf("%.3f", alpha_2)))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16)
    ) +
    ylim(1, 9)
  
  c(alpha_lambda1 = alpha_1, alpha_lambda2 = alpha_2)
}

alpha_diag <- estimate_alpha_diag(H1, H2, h_values)
alpha_eigen <- estimate_alpha_eigen(H1, H2, h_values)

subtitle_diag <- textGrob("Diagonal Elements", x = 0.418, hjust = 0, gp = gpar(fontsize = 16))
subtitle_eigen <- textGrob("Eigenvalues", x = 0.455, hjust = 0, gp = gpar(fontsize = 16))

grid.arrange(
  subtitle_diag,
  arrangeGrob(p1, p2, ncol = 2),
  textGrob(""),
  subtitle_eigen,
  arrangeGrob(p4, p3, ncol = 2),
  ncol = 1,
  heights = c(0.1, 1, 0.03, 0.1, 1)
)
