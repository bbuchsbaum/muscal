library(testthat)
library(muscal)

context("BaMFA reconstruction")

test_that("BaMFA reconstructs blocks with low error", {
  set.seed(1)
  Xs <- lapply(1:4, function(i) matrix(rnorm(40), 10, 4) + i)
  names(Xs) <- paste0("Block_", 1:4)
  
  fit <- bamfa(Xs, k_g = 2L, k_l = 1L, niter = 5L)
  recon <- predict(fit, Xs)
  
  # check that predict with new_data is reasonably close to without
  recon_train <- predict(fit)
  expect_equal(recon, recon_train, tolerance=0.05)
  
  
  err   <- mean(vapply(seq_along(Xs),
                       function(i) {
                         mean((Xs[[i]] - recon[[i]])^2)
                       },
                       numeric(1)))
                       
  # A basic sanity check that error is not enormous
  expect_lt(err, 1.5)
})

test_that("BaMFA reconstruction error is low for a known structure", {
  set.seed(42)
  n_obs <- 20
  p1 <- 15
  p2 <- 20
  p3 <- 10
  
  # Global structure
  G_true <- qr.Q(qr(matrix(rnorm( (p1+p2+p3) * 2), p1+p2+p3, 2)))
  S_true <- matrix(rnorm(n_obs * 2), n_obs, 2)
  
  # Local structure
  B1_true <- qr.Q(qr(matrix(rnorm(p1*1), p1, 1)))
  B2_true <- qr.Q(qr(matrix(rnorm(p2*1), p2, 1)))
  B3_true <- qr.Q(qr(matrix(rnorm(p3*1), p3, 1)))
  U1_true <- matrix(rnorm(n_obs*1), n_obs, 1)
  U2_true <- matrix(rnorm(n_obs*1), n_obs, 1)
  U3_true <- matrix(rnorm(n_obs*1), n_obs, 1)
  
  X1 <- S_true %*% t(G_true[1:p1,]) + U1_true %*% t(B1_true) + matrix(rnorm(n_obs*p1, sd=0.1), n_obs, p1)
  X2 <- S_true %*% t(G_true[(p1+1):(p1+p2),]) + U2_true %*% t(B2_true) + matrix(rnorm(n_obs*p2, sd=0.1), n_obs, p2)
  X3 <- S_true %*% t(G_true[(p1+p2+1):(p1+p2+p3),]) + U3_true %*% t(B3_true) + matrix(rnorm(n_obs*p3, sd=0.1), n_obs, p3)
  
  Xs <- list(X1=X1, X2=X2, X3=X3)
  
  fit <- bamfa(Xs, k_g = 2L, k_l = 1L, niter = 15L, tol=1e-7)
  recon <- predict(fit)
  
  err   <- mean(vapply(seq_along(Xs),
                       function(i) {
                         sum((Xs[[i]] - recon[[i]])^2) / sum(Xs[[i]]^2)
                       },
                       numeric(1)))
                       
  # Expect relative error to be reasonably small
  expect_lt(err, 0.25)
}) 