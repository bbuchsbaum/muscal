library(testthat)
library(muscal)
library(multidesign)

skip_if_not_installed("multidesign")

test_that("bada.multidesign works with basic splitting", {
  set.seed(123)
  x <- matrix(rnorm(20), nrow = 10)
  design <- data.frame(
    subj_id = rep(c("S1", "S2"), each = 5),
    y = rep(c("A", "B"), 5)
  )
  md <- multidesign::multidesign(x, design)
  
  res <- bada(md, y = y, subject = subj_id, ncomp = 1)
  
  expect_equal(levels(res$subjects), c("S1", "S2"))
})

