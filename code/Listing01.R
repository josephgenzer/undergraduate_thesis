n <- 10000
mu <- 5
lambda <- 1 / mu
sigma <- 2

normal_samples <- rnorm(n, mean = mu, sd = sigma)
exp_samples <- rexp(n, rate = lambda)

cummean_normal <- cumsum(normal_samples) / seq_len(n)
cummean_exp <- cumsum(exp_samples) / seq_len(n)

plot(cummean_normal, type = "l", col = "#4682B4", lwd = 1.5,
     ylim = c(mu - 1, mu + 1),
     xlab = "Sample Size (n)",
     ylab = expression("Sample Mean (" * bar(X)[n] * ")"))
lines(cummean_exp, col = "#B22222", lwd = 1.5)
legend("topright",
       legend = c("Normal", "Exponential"),
       col = c("#4682B4", "#B22222"),
       lwd = 1.5,
       bg = "white",
       box.lwd = 1)
