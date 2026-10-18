test_that("loading reliability matches signs and component permutations", {
  reference <- diag(4)[, 1:3]
  draw <- sweep(reference[, c(3, 1, 2)], 2, c(-1, 1, -1), "*")
  result <- loading_reliability(list(sdev = c(3, 2, 1)), V_ref = reference,
                                 boot_loadings = list(reference, draw))
  expect_equal(result$mean, reference)
  expect_equal(result$sd, matrix(0, 4, 3))
  expect_equal(result$lower, reference)
  expect_equal(result$upper, reference)
  expect_equal(result$alignment$permutations[[2]], c(2L, 3L, 1L))
})

test_that("assignment maximizes total matching rather than greedy columns", {
  cost <- matrix(c(1, 2, 100, 2, 100, 100, 100, 3, 1), 3, 3, byrow = TRUE)
  assignment <- muscal:::.muscal_min_assignment(cost)
  expect_equal(assignment, c(2L, 1L, 3L))
  expect_equal(sum(cost[cbind(1:3, assignment)]), 5)
  set.seed(342)
  permutations <- as.matrix(expand.grid(rep(list(1:4), 4)))
  permutations <- permutations[apply(permutations, 1, function(p) length(unique(p)) == 4), ]
  for (iteration in 1:4) {
    cost <- matrix(runif(16, -2, 5), 4, 4)
    actual <- muscal:::.muscal_min_assignment(cost)
    enumerated <- apply(permutations, 1, function(p) sum(cost[cbind(1:4, p)]))
    expect_equal(sum(cost[cbind(1:4, actual)]), min(enumerated), tolerance = 1e-12)
    expect_equal(sort(actual), 1:4)
  }
})

test_that("tied-root rotations have subspace evidence but unavailable coordinate intervals", {
  reference <- diag(4)[, 1:3]
  rotate <- function(angle) {
    R <- diag(3)
    R[1:2, 1:2] <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2, 2)
    reference %*% R
  }
  fit <- list(sdev = c(3, 3 * (1 - 1e-9), 1))
  result <- loading_reliability(fit, V_ref = reference,
    boot_loadings = list(rotate(pi / 6), rotate(-pi / 9)))
  for (field in c("mean", "sd", "lower", "upper", "sign_consistency")) {
    expect_true(all(is.na(result[[field]][, 1:2])))
  }
  expect_equal(result$mean[, 3], reference[, 3])
  expect_equal(result$alignment$ambiguous_components, c(TRUE, TRUE, FALSE))
  expect_equal(as.numeric(result$alignment$subspace_distance), c(0, 0), tolerance = 1e-12)
})

test_that("loading matching rejects invalid bases and reports unknown identifiability", {
  reference <- diag(3)[, 1:2]
  unknown <- loading_reliability(list(), V_ref = reference,
                                boot_loadings = list(reference, reference))
  expect_identical(unknown$alignment$identifiability, "unknown_reference_spectrum")
  bad <- reference
  bad[, 1] <- 0
  expect_error(loading_reliability(list(sdev = c(3, 1)), V_ref = reference,
    boot_loadings = list(bad, reference)), "zero.norm")
  bad[1, 1] <- Inf
  expect_error(loading_reliability(list(sdev = c(3, 1)), V_ref = reference,
    boot_loadings = list(bad, reference)), "finite")
  expect_error(loading_reliability(list(sdev = diag(2)), V_ref = reference,
    boot_loadings = list(reference, reference)), "spectrum.*vector")
})

test_that("loading alignment matches named feature rows without guessing identities", {
  reference <- diag(3)[, 1:2]
  rownames(reference) <- paste0("f", 1:3)
  permuted <- reference[3:1, ]
  result <- loading_reliability(list(sdev = c(3, 1)), V_ref = reference,
                                 boot_loadings = list(permuted, reference))
  expect_equal(result$mean, unname(reference))
  rownames(permuted) <- paste0("other", 1:3)
  expect_error(loading_reliability(list(sdev = c(3, 1)), V_ref = reference,
    boot_loadings = list(permuted, reference)), "feature identities")
})
