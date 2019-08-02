test_that("normalization works as expected", {
  n_prots <- 1000
  mu_g <- rnorm(n_prots, mean=25, sd=3)
  samp1 <- rnorm(n_prots, mu_g+0.4, sd=0.3)
  samp2 <- rnorm(n_prots, mu_g-0.1, sd=0.3)
  samp3 <- rnorm(n_prots, mu_g+0.2, sd=0.3)
  samp4 <- rnorm(n_prots, mu_g-0.6, sd=0.3)
  X <- cbind(samp1, samp2, samp3, samp4)
  
  # plot(X[, 1], X[, 1] - rowMeans(X)); abline(h = median(X[, 1] - rowMeans(X)))
  
  Xnorm <- median_normalization(X)
  
  # plot(Xnorm[, 1], Xnorm[, 1] - rowMeans(Xnorm)); abline(h = median(Xnorm[, 1] - rowMeans(Xnorm)))
  
  new_meds <- vapply(1:4, function(idx) median(Xnorm[, idx] - rowMeans(Xnorm)), FUN.VALUE = 0.0)
  expect_true(all(abs(new_meds) < 0.01))
})





test_that("normalization works with consistent bias using spike in correction ", {
  n_spike_in <- 50
  n_prots <- 1000
  n_off <- 250
  mu_g <- rnorm(n_prots + n_off + n_spike_in, mean=25, sd=3)
  samp1 <- c(rnorm(n_prots, mu_g[seq_len(n_prots)]+0.4, sd=0.3), rnorm(n_off, mu_g[n_prots + seq_len(n_off)] + 1, sd=0.3), rnorm(n_spike_in, mu_g[n_prots + n_off + seq_len(n_spike_in)]+0.4, sd=0.01))
  samp2 <- c(rnorm(n_prots, mu_g[seq_len(n_prots)]-0.1, sd=0.3), rnorm(n_off, mu_g[n_prots + seq_len(n_off)] + 1, sd=0.3), rnorm(n_spike_in, mu_g[n_prots + n_off + seq_len(n_spike_in)]-0.1, sd=0.01))
  samp3 <- c(rnorm(n_prots, mu_g[seq_len(n_prots)]+0.2, sd=0.3), rnorm(n_off, mu_g[n_prots + seq_len(n_off)] - 1, sd=0.3), rnorm(n_spike_in, mu_g[n_prots + n_off + seq_len(n_spike_in)]+0.2, sd=0.01))
  samp4 <- c(rnorm(n_prots, mu_g[seq_len(n_prots)]-0.6, sd=0.3), rnorm(n_off, mu_g[n_prots + seq_len(n_off)] - 1, sd=0.3), rnorm(n_spike_in, mu_g[n_prots + n_off + seq_len(n_spike_in)]-0.6, sd=0.01))
  X <- cbind(samp1, samp2, samp3, samp4)
  
  # plot(X[, 1], X[, 1] - rowMeans(X), col=ifelse(seq_along(mu_g) <= n_prots, "black", "red"),
  #      pch = 16, cex=0.8)
  # abline(h = median(X[seq_len(n_prots),  1] - rowMeans(X[seq_len(n_prots), ])))
  
  Xnorm <- median_normalization(X, n_prots + n_off + seq_len(n_spike_in))
  
  # plot(Xnorm[, 1], Xnorm[, 1] - rowMeans(Xnorm), col=ifelse(seq_along(mu_g) <= n_prots, "black", "red"),
  #      pch = 16, cex=0.8)
  # abline(h = median(Xnorm[seq_len(n_prots),  1] - rowMeans(Xnorm[seq_len(n_prots), ])))
  
  new_meds <- vapply(1:4, function(idx) median(Xnorm[seq_len(n_prots), idx] - rowMeans(Xnorm[seq_len(n_prots), ])), FUN.VALUE = 0.0)
  expect_true(all(abs(new_meds) < 0.1))
})
